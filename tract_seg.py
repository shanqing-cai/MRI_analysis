#!/usr/bin/python

import os
import sys
import glob
import time
import argparse
from subprocess import Popen, PIPE
import numpy as np

import nibabel as nib

from aparc12_probtrackx_2 import get_roi_num
from aparc12 import get_aparc12_cort_rois
from scai_utils import *

# import slaparc

# ==== Definition: get_region_mean ==== #
def get_region_mean(roiNames, roiSizes, \
                    maskedMeans, regions, regionROIs, region):
    if regions.count(region) == 0:
        raise IOError, 'Region %s not found.'%region
    idx = regions.index(region)
    inc_rois = regionROIs[idx]

    t_roiSizes = []
    t_maskedMeans = []
    for i0 in range(len(roiNames)):
        if list(inc_rois).count(roiNames[i0]) == 1:
            t_roiSizes.append(roiSizes[i0])
            t_maskedMeans.append(maskedMeans[i0])

    #print(t_roiSizes)
    #print(t_maskedMeans)

    denom = np.sum(np.array(t_roiSizes))
    if denom > 0:
        region_total = np.sum(np.array(t_roiSizes) * np.array(t_maskedMeans)) \
                       / denom;
    else:
        region_total = 0
    return region_total
# ==== ~Definition: get_region_mean ==== #

# Config: directories:
tractSegDir = '/users/cais/STUT/analysis/tractseg_aparc12/'
traculaDir = '/users/cais/STUT/analysis/dti2/tracula'
FSDATADir = '/users/cais/STUT/FSDATA'

ROI_MASK_BASE = "/users/cais/STUT/analysis/aparc12_tracts_2"

CORTICAL_CTAB = "/users/cais/STUT/slFRS17.ctab"
SUBCORTICAL_CTAB = "/software/atlas/ASAP_subcortical_labels.txt"

ALL_STEPS = ["track", "stats", "seg", "zip"]

# ==== CONFIG: Construct the parcellation profile ==== #
leaveOutROIs = ['pINS', 'aINS', 'LG']

regions = ['Prefrontal', 'Premotor', 'Precentral', 'Postcentral', \
           'PPC', 'Occipital', 'Temporal']
regionCodes = [1, 2, 3, 4, 5, 6, 7]
"""
regionROIs = [['SFg', 'aMFg', 'iFo', 'iFt', 'FMC', 'FO', 'FOC', 'FP', \
               'pMFg'], \
              ['pSMA', 'aSMA', 'adPMC', 'mdPMC', 'pdPMC', 'vPMC'], \
              ['aCO', 'vMC', 'dMC'], \
              ['vSSC', 'pCO', 'dSSC'], \
              ['AG', 'aSMg', 'pSMg', 'SPL'], \
              ['CALC', 'MTO', 'ITO'], \
              ['H', 'PT', 'aSTg', 'aMTg', 'aITg', 'pSTg', 'pMTg', 'pITg', \
               'TP', 'adSTs', 'avSTs', 'pdSTs', 'pvSTs'],\
             ]
"""
regionROIs = [[]] * len(regions)
for (i0, region) in enumerate(regions):
    regionROIs[i0] = get_aparc12_cort_rois(lobe=region, bSpeech=False)

allROIs = get_aparc12_cort_rois(lobe="all", bSpeech=False)

# ==== ~CONFIG: Construct the parcellation profile ==== #


# ==== sub_routine ==== #
def tract_seg_sub(args, step, roiName=""):
    # Check input arguments
    sID = args.sID
    seedROI = args.seedROI
    b_ccStop = args.b_ccStop

    if ALL_STEPS.count(step) == 0 and step != "roiconn":
        raise Exception, "Unrecognized step: %s" % step

    if step == "roiconn":
        if allROIs.count(roiName) == 0:
            raise Exception, \
                  "Under step == %s, found unrecognized ROI name %s" \
                  % (step, roiName)
    
    # If the ccMask is needed, locate or generate the ccMask:
    if b_ccStop:
        ccMask_fn = os.path.join(FSDATADir, sID, 'mri', \
                                 'ccMask_diff_aparc12.nii.gz')
        if not os.path.isfile(ccMask_fn):
            saydo('python /users/cais/STUT/scripts/gen_ccMask_aparc12.py %s' \
                  % sID)
            if not os.path.isfile(ccMask_fn):
                raise IOError, 'Failed to generate ccMask.'
        else:
            print('\nINFO: ccMask_fn = %s' % ccMask_fn)


    # Locate the diff-space ROI file
    if len(seedROI) <= 3:
        raise IOError, 'Unrecognized seed ROI: %s'%seedROI

    if seedROI[:3] != 'lh.' and seedROI[:3] != 'rh.':
        raise IOError, 'Unrecognized seed ROI: %s'%seedROI

    hemi = seedROI[:2]

    # === Generate diffusion-space seed mask === %
    check_file(SUBCORTICAL_CTAB)
    ctab = read_ctab(SUBCORTICAL_CTAB)
    seedROI1 = seedROI.replace("lh.", "Left-").replace("rh.", "Right-")

    if ctab[1].count(seedROI1) != 1:
        raise Exception, "Cannot find entry for ROI %s in color table file %s" \
                         % (seedROI, SUBCORTICAL_CTAB)
    roiCode = ctab[0][ctab[1].index(seedROI1)]

    print("\nINFO: roiCode for %s = %d" % (seedROI, roiCode))

    diffROI = os.path.join(FSDATADir, sID, 'label', 'aparcSL_' + hemi, \
                           seedROI + '.diff_aparc12.nii.gz')
    
    if not os.path.isfile(diffROI):
        # TODO: use gen_subCortMask_aparc12.py
        gen_roiMask_cmd = "gen_subCortMask_aparc12.py %s %d %s" \
                          % (sID, roiCode, diffROI)
        saydo(gen_roiMask_cmd)
        check_file(diffROI)

    print("\nINFO: dffusion-space seed mask file: %s" % diffROI)

    # === Generate ASCII-format coordinates file === #
    if b_ccStop:
        roiDir = os.path.join(tractSegDir, sID, seedROI + '_ccStop')
    else:
        roiDir = os.path.join(tractSegDir, sID, seedROI)

    check_dir(roiDir, bCreate=True)

    coordFN = os.path.join(roiDir, '%s_coords.txt' % seedROI)
    if not os.path.isfile(coordFN):        
        delete_file_if_exists(coordFN)
        matlabCmd = 'matlab -nodesktop -nosplash -r "img_nz_coord ' \
                    + diffROI + \
                    ' ' + coordFN + '; exit"'
        saydo(matlabCmd)

        check_file(coordFN)

    # === Read all voxel coordinates (subs) in the ROI === # 
    coord_f = open(coordFN, "rt")
    coord_t = coord_f.read()
    coord_f.close()

    coord_t = coord_t.split('\n')
    
    coords = []


    for t_lin in coord_t:
        if len(t_lin) == 0:
            continue;
        
        idx0 = t_lin.index('[')
        idx1 = t_lin.index(']')
        t_txt = t_lin[(idx0 + 1) : idx1]

        if t_txt.count(', ') == 2:
            coords.append([int(t_txt.split(', ')[0]), \
                           int(t_txt.split(', ')[1]), \
                           int(t_txt.split(', ')[2])])

    if step == "track":
        # Create all the mask files and the command
        probtrackxCmds = []
        jobNames = []
        expectFNs = []
        seedFNs = []

        merged = os.path.join(traculaDir, sID, 'dmri.bedpostX', 'merged')

        mask = os.path.join(traculaDir, sID, 'dmri', 'nodif_brain_mask.nii.gz')
        check_file(mask)

        check_dir(os.path.join(roiDir, 'vox'), bCreate=True)

        for coord in coords:
            outputDir = os.path.join(roiDir, '%d_%d_%d'%(coord[0], coord[1], coord[2]))
            t_expectFN = os.path.join(outputDir, 'fdt_paths.nii.gz');
    
            seedFN = os.path.join(roiDir, 'vox', 'vox_%d_%d_%d.diff.nii.gz'\
                                  %(coord[0], coord[1], coord[2]))
        
            if not os.path.isfile(t_expectFN):            
                if not os.path.isfile(seedFN):
                    gen_mask_cmd = 'fslmaths %s -roi %d 1 %d 1 %d 1 0 1 %s'\
                           %(diffROI, coord[0], coord[1], coord[2], seedFN)
    
                    saydo(gen_mask_cmd)
                    check_file(seedFN)
            
                    seedFNs.append(seedFN)

                if b_ccStop:
                    # probtrackxCmds.append('probtrackx2 --mode=seedmask -x %s -s %s -m %s -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd --pd --stop=%s --dir=%s' % (seedFN, merged, mask, ccMask_fn, outputDir))
                    probtrackxCmds.append('probtrackx2 -x %s -s %s -m %s -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd --pd --stop=%s --dir=%s' % (seedFN, merged, mask, ccMask_fn, outputDir))

                else:
                    # probtrackxCmds.append('probtrackx2 --mode=seedmask -x %s -s %s -m %s -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd --dir=%s'%(seedFN, merged, mask, outputDir))
                    probtrackxCmds.append('probtrackx2 -x %s -s %s -m %s -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd --dir=%s'%(seedFN, merged, mask, outputDir))

                jobNames.append('%d_%d_%d'%(coord[0], coord[1], coord[2]))
                expectFNs.append(t_expectFN)
    
        for (i0, cmd) in enumerate(probtrackxCmds):
            saydo(cmd)
            check_file(expectFNs[i0])
    elif step == "stats":
        # === Load the maskROIs === #
        maskROIs_fn = glob.glob(os.path.join(ROI_MASK_BASE, sID, \
                         "aparc12_%s_*.diff.nii.gz" % hemi))
        maskROIs_fn.sort()

        maskROIs = []        

        maskROINums = []
        for fn in maskROIs_fn:
            [foo, t_fn] = os.path.split(fn)
            
            t_ROI_name = t_fn.split(".")[0].split("_")[-1]

            maskROIs.append(t_ROI_name)
            maskROINums.append(get_roi_num(CORTICAL_CTAB, hemi, t_ROI_name))

        # === Get the sizes of the ROIs === #
        nROIs = len(maskROIs)
        roiNVoxes = [float('nan')] * nROIs
        for i0 in range(len(maskROIs_fn)):
            t_mask = maskROIs_fn[i0]
            (stdout, stderr) = Popen(['fslstats', t_mask, '-V'], stdout=PIPE)\
                .communicate()
            roiNVoxes[i0] = float(stdout.split(' ')[0])

        # === Process the voxels, one by one === %
        for (k0, coord) in enumerate(coords):
            print("Step %s: processing voxel #%d of %d..." \
                  % (step, k0, len(coords)))
            
            voxDir = os.path.join(roiDir, \
                                  '%d_%d_%d'%(coord[0], coord[1], coord[2]))
            check_dir(voxDir)
            masked_mean_fn = os.path.join(voxDir, 'masked_mean.txt')
            max_roi_fn = os.path.join(voxDir, 'max_roi.txt')

            if True:
            #if (not os.path.isfile(masked_mean_fn)) or \
            #   (not os.path.isfile(max_roi_fn)):
                maskedVals = [float('nan')] * nROIs

                fdtPathsFN = os.path.join(voxDir, 'fdt_paths.nii.gz')
                check_file(fdtPathsFN)

                if not os.path.isfile(fdtPathsFN):
                    print("WARNING: fdt_paths.nii.gz not found, skipping directory %s" % voxDir)
                    continue

                t_txt = ''
                for i0 in range(len(maskedVals)):
                    t_mask = maskROIs_fn[i0]

                    #(stdout, stderr) = Popen(['fslstats', t_mask, '-V'], \
                    #                         stdout=PIPE, stderr=PIPE)\
                    #                  .communicate()
                    #maskNVox = int(stdout.split(" ")[0])                    

                    if roiNVoxes[i0] == 0:
                        maskedVals[i0] = 0.0;
                        print("WARNING: ROI mask %s is empty." % t_mask)
                    else:                        
                        (stdout, stderr) = Popen(['fslstats', fdtPathsFN, \
                                              '-k', t_mask, \
                                              '-m'], stdout=PIPE, stderr=PIPE)\
                                              .communicate()

                        if len(stderr) > 0:
                            raise Exception, "fslstats reported ERROR during processing of coordinate (%d, %d, %d) and cortical ROI %s: \n%s" \
                               % (coord[0], coord[1], coord[2], t_mask, stderr)

                        maskedVals[i0] = float(stdout.split(' ')[0])
                        
                    #print("roiNVoxes = %d; maskNVox = %d" \
                    #      % (roiNVoxes[i0], maskNVox)) # DEBUG
                    t_txt += '%s %f %f\n' \
                         % (maskROIs[i0], roiNVoxes[i0], maskedVals[i0])        

                masked_mean_f = open(masked_mean_fn, "wt")
                masked_mean_f.write(t_txt)
                masked_mean_f.close()
                print('Wrote masked means to %s'%masked_mean_fn)

                maxVal = np.max(maskedVals)
                maxIdx = maskedVals.index(maxVal)
                maxROI = maskROIs[maxIdx]
    
                max_roi_f = open(max_roi_fn, "wt")
                max_roi_f.write('%s %d %f %f\n'\
                                %(maxROI, maskROINums[maxIdx], \
                                roiNVoxes[maxIdx], maxVal))
                max_roi_f.close()
                print('Wrote max ROI to %s\n'%max_roi_fn)

    elif step == "seg":
        #sys.exit(0)
        add_cmd = 'fslmaths '
        add_reg_cmd = 'fslmaths '

        a_region_means = np.zeros([len(coords), len(regions)])
        add_regMean_cmd = ["fslmaths "] * len(regions)
        
        #print(add_regMean_cmd[6])
        #sys.exit(0)

        for (j0, coord) in enumerate(coords):
            voxDir = os.path.join(roiDir, \
                                  '%d_%d_%d'%(coord[0], coord[1], coord[2]))
            check_dir(voxDir)
            
            masked_mean_fn = os.path.join(voxDir, 'masked_mean.txt')
            check_file(masked_mean_fn)

            max_roi_fn = os.path.join(voxDir, 'max_roi.txt')
            check_file(max_roi_fn)
            
            masked_mean_f = open(masked_mean_fn, 'r')
            masked_mean_t = masked_mean_f.read()
            masked_mean_f.close()            

            mean_t = masked_mean_t.split('\n');
            mean_t = remove_empty_strings(mean_t)

            if len(mean_t) > len(allROIs):
                raise Exception, "Number of lines in file %s (%d) is greater than the total number of aparc12 ROIs (%d)" \
                                 % (masked_mean_fn, len(mean_t), len(allROIs))

            t_roiNames = []
            t_roiSizes = []
            t_roiNums = []
            t_maskedMeans = []
            
            for t_lin in mean_t:
                if len(t_lin) == 0:
                    continue
                t_items = t_lin.split(' ')
                t_roiNames.append(t_items[0])
                t_roiSizes.append(float(t_items[1]))
                t_roiNums.append(get_roi_num(CORTICAL_CTAB, hemi, t_items[0]))
                t_maskedMeans.append(float(t_items[2]))
            
            for removeItem in leaveOutROIs:
                if t_roiNames.count(removeItem) == 1:
                    idx = t_roiNames.index(removeItem)
                    t_roiNames.pop(idx)
                    t_roiNums.pop(idx)
                    t_maskedMeans.pop(idx)

            #print(t_maskedMeans)
            #sys.exit(0)
    
            # ROI labeling
            max_roi_num = t_roiNums[t_maskedMeans.index(np.max(t_maskedMeans))]
            max_roi_name = \
                t_roiNames[t_maskedMeans.index(np.max(t_maskedMeans))]
    
            print('Voxel (%d, %d, %d): max ROI = %s, num = %d'\
                      %(coord[0], coord[1], coord[2], \
                        max_roi_name, max_roi_num))

            voxImgFN = os.path.join(roiDir, "vox", \
                                    'vox_%d_%d_%d.diff.nii.gz'\
                                    % (coord[0], coord[1], coord[2]))
            voxImgFN_mul = os.path.join(roiDir, "vox", \
                                        'vox_%d_%d_%d_mul.diff.nii.gz'\
                                        % (coord[0], coord[1], coord[2]))

            mult_cmd = 'fslmaths %s -mul %d %s' \
                       % (voxImgFN, max_roi_num, voxImgFN_mul)
            delete_file_if_exists(voxImgFN_mul)
            saydo(mult_cmd)
            check_file(voxImgFN_mul)

            if add_cmd.count(' ') == 1:
                add_cmd += voxImgFN_mul + ' '
            else:
                add_cmd += '-add ' + voxImgFN_mul + ' '
    
            # Region labels
            t_region_means = []
            for region in regions:
                t_region_means.append(get_region_mean(t_roiNames, t_roiSizes, \
                              t_maskedMeans, regions, regionROIs, region))

            max_region_num = \
                regionCodes[t_region_means.index(np.max(t_region_means))]
            max_region_name = \
                regions[t_region_means.index(np.max(t_region_means))]
    
            print('\tmax region = %s, num = %d' \
                  % (max_region_name, max_region_num))         
            
            #print(t_region_means)
            #print(a_region_means)
            #print(add_cmd)
            #sys.exit(0)
    
            voxImgFN_reg = os.path.join(roiDir, \
                                        'vox_%d_%d_%d_reg.diff.nii.gz'\
                                        %(coord[0], coord[1], coord[2]))
            mult_reg_cmd = 'fslmaths %s -mul %d %s'\
                           %(voxImgFN, max_region_num, voxImgFN_reg)

            delete_file_if_exists(voxImgFN_reg)
            saydo(mult_reg_cmd)
            check_file(voxImgFN_reg)

            if add_reg_cmd.count(' ') == 1:
                add_reg_cmd += voxImgFN_reg + ' '
            else:
                add_reg_cmd += '-add ' + voxImgFN_reg + ' '

            # Regional means
            for k1 in range(len(regions)):
                a_region_means[j0][k1] = t_region_means[k1]
                vox_regMean_fn = os.path.join(roiDir, \
                                 "vox_%d_%d_%d_regMean_%d.diff.nii.gz" \
                                 % (coord[0], coord[1], coord[2], k1))
                #print(vox_regMean_fn)

                mult_regMean_cmd = "fslmaths %s -mul %f %s" \
                           % (voxImgFN, t_region_means[k1], vox_regMean_fn)
                saydo(mult_regMean_cmd)
                check_file(vox_regMean_fn)

                #if add_regMean_cmd[k1].count(' ') == 1:
                #    add_regMean_cmd[k1] += vox_regMean_fn + ' '
                #else:
                #    add_regMean_cmd[k1] += '-add ' + vox_regMean_fn + ' '
                
            #if j0 == 1:
            #    print(add_regMean_cmd[0])
            #    print(add_regMean_cmd[6])
            #    sys.exit(0)

        tract_seg_fn = os.path.join(roiDir, 'tract_parc.nii.gz')
        tract_region_seg_fn = os.path.join(roiDir, 'tract_regionParc.nii.gz')

        #add_cmd += tract_seg_fn
        #add_reg_cmd += tract_region_seg_fn

        #delete_file_if_exists(tract_seg_fn)
        #saydo(add_cmd)
        #check_file(tract_seg_fn)

        add_cmd = 'glob_fsl_add.py "%s" %s -v --cleanup' \
                  % (os.path.join(roiDir, "vox_*_*_*_reg.diff.nii.gz"), \
                     tract_region_seg_fn)

        saydo(add_cmd)
        check_file(tract_region_seg_fn)

        print("\nINFO: lobe-by-lobe segmentation file saved at\n\t%s" \
              % tract_region_seg_fn)

        for k1 in range(len(regions)):
            tract_region_mean_fn = \
                os.path.join(roiDir, "tract_regionMean_%d.nii.gz" % k1)

            add_cmd_regMean = 'glob_fsl_add.py "%s" %s -v --cleanup' \
                  % (os.path.join(roiDir, \
                                  "vox_*_*_*_regMean_%d.diff.nii.gz" % k1), \
                     tract_region_mean_fn)
            saydo(add_cmd_regMean)
            check_file(tract_region_mean_fn)
            
            print("\nINFO: lobe-by-lobe mean tract density saved at\n\t%s" \
                  % tract_region_mean_fn)
            

        #delete_file_if_exists(tract_region_seg_fn)
        #saydo(add_reg_cmd)
        #check_file(tract_region_seg_fn)

        #print("\nINFO: ROI-by-ROI segmentation file saved at\n\t%s" \
        #      % tract_seg_fn)

    elif step == "roiconn":
        print("INFO: step = roiconn; roiName = %s" % roiName)
        
        voxMaskDir = os.path.join(roiDir, "vox")
        check_dir(voxMaskDir)

        voxMask0 = os.path.join(voxMaskDir, \
                                "vox_%d_%d_%d.diff.nii.gz" \
                                % (coords[0][0], coords[0][1], coords[0][2]))
        check_file(voxMask0)

        img = nib.load(voxMask0)
        imgdat = img.get_data()

        for (k0, coord) in enumerate(coords):
            #print("Step %s: processing voxel #%d of %d..." \
            #      % (step, k0, len(coords)))
            
            voxDir = os.path.join(roiDir, \
                                  '%d_%d_%d'%(coord[0], coord[1], coord[2]))
            check_dir(voxDir)
            masked_mean_fn = os.path.join(voxDir, 'masked_mean.txt')

            mmtxt = read_text_file(masked_mean_fn)
            # == Locate the ROI line and get the masked mean == #
            bROIFound = False
            for (j0, t_line) in enumerate(mmtxt):
                if len(mmtxt) == 0:
                    continue
                t_items = t_line.split(" ")
                if len(t_items) != 3:
                    continue
                
                if t_items[0] == roiName:
                    bROIFound = True
                    break
            
            if not bROIFound:
                raise Exception, "Failed to find entry for ROI %s in file %s" \
                      % (roiName, masked_mean_fn)
            t_maskedMean = float(t_items[-1])

            imgdat[coord[0], coord[1], coord[2]] = t_maskedMean
            #print("Processing voxel %d/%d: Assigning value [%d, %d, %d] --> %s" \
            #      % (k0, len(coords), coord[0], coord[1], coord[2], \
            #         imgdat[coord[0], coord[1], coord[2]]))

        # Write to output file: 
        roiconn_dir = os.path.join(roiDir, "aparc12")
        check_dir(roiconn_dir, bCreate=True)

        roiConnFN = os.path.join(roiconn_dir, "%s.nii.gz" % roiName)
        imgout = nib.nifti1.Nifti1Image(imgdat, img.get_affine(), \
                                        header=img.get_header())
        imgout.to_filename(roiConnFN)
        check_file(roiConnFN)

        print("INFO: ROI connectivity file saved at: %s" % (roiConnFN))
        
    elif step == "zip":
        tract_region_seg_fn = os.path.join(roiDir, 'tract_regionParc.nii.gz')
        
        check_file(tract_region_seg_fn)
        dvds0 = glob.glob(os.path.join(roiDir, "??_??_??"))
        dvds = []
        for t_dvd in dvds0:
            if os.path.isdir(t_dvd):
                dvds.append(t_dvd)

        tarOut = os.path.join(roiDir, "intermediate.tar.gz")        
        if len(dvds) == 0:
            if os.path.isfile(tarOut):
                print("It appears that the zip step has already been performed")
                return
            else:
                raise Exception, "Unexpected content in directory: %s" % roiDir
        
        print("INFO: Found %d directories to zip" % (len(dvds) + 1))
        
        cwd = os.getcwd()
        os.chdir(roiDir)
        tarCmd = "tar czvf %s %s " % (tarOut, "vox")
        for dvd in dvds:
            (foo, dn) = os.path.split(dvd)
            tarCmd += "%s " % dn
        saydo(tarCmd)

        rmCmd = tarCmd.replace("tar czvf %s" % tarOut, "rm -rf")
        saydo(rmCmd)

        os.chdir(cwd)
    
    return
# ==== ~sub_routine ==== #


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Parcellation of the subcortical ROIs based on probabilistic tractograpy\nAuthor: Shanqing Cai (shanqing.cai@gmail.com) \nDate: 2013-02-27")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("seedROI", help="seedROI")
    ap.add_argument("step", help="Step {track, stats, seg, zip, roiconn or all {all=tract+stats+seg+zip}}")
    ap.add_argument("--ccStop", dest="b_ccStop", action="store_true", \
                    help="Use corpus callosum (CC) as a stop mask")
    ap.add_argument("--roi", dest="roi", type=str, default="", \
                    help="ROI name (for use with step==roiconn")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    # === Parse input argument === #
    args = ap.parse_args()
    step = args.step
    
    if step == "roiconn":
        roi = args.roi
        assert(roi != "")

    # === === #
    print("INFO: step = %s" % step)
    if step == "all":
        for (i0, t_step) in enumerate(ALL_STEPS):
            tract_seg_sub(args, t_step)
    else:
        if step != "roiconn":
            tract_seg_sub(args, step)
        else:
            tract_seg_sub(args, step, roiName=roi)
        
