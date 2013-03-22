#!/usr/bin/python

import os
import sys
import glob
import pickle
import argparse
import numpy as np
import scipy.stats as stats
from subprocess import Popen, PIPE
from scai_utils import read_ctab
#import matplotlib.pyplot as plt

DATA_dir = '/users/cais/STUT/DATA'
FSDATA_dir = '/users/cais/STUT/FSDATA'
ANALYSIS_DIR = '/users/cais/STUT/analysis'
bips_resting_dir = '/users/cais/STUT/analysis'
bips_resting_dir_2 = "/users/cais/STUT/analysis/resting_bips_2"
#ASAP_TABLE = '/software/atlas/ASAP_labels.txt'
APARC12_TABLE = '/users/cais/STUT/slFRS17.ctab'
ASAP_SC_TABLE = '/software/atlas/ASAP_subcortical_labels.txt'

HEMIS = ['lh', 'rh']

def saydo(cmd):
    print('\n%s\n'%cmd)
    os.system(cmd)


def get_roi_ids(asap_table, sc_table):
    #=== Cortical ROIs ===#
    tablef = open(asap_table, 'r')
    txt = tablef.read()
    
    tablef.close()
    txt = txt.split('\n')

    t_rois = []
    t_ids = []
    for t in txt:
        if len(t) == 0:
            continue

        tt = t.split(' ')

        while tt.count('') > 0:
            tt.remove('')

        if len(tt) == 0:
            continue

        if tt[1] == 'None' or tt[1] == 'White' or tt[1] == 'Gray' \
                or tt[1] == 'CN' or tt[1].startswith('None') \
                or tt[1] == 'Unknown':
            continue

        t_rois.append(tt[1])
        t_ids.append(tt[0])
    
    s_rois = t_rois
    s_ids = t_ids

    b_rois = []
    b_ids = []

    #=== Subcortical ROIs ===#
    tablef = open(sc_table, 'r')
    txt = tablef.read()
    tablef.close()
    txt = txt.split('\n')

    t_rois = []
    t_ids = []
    for t in txt:
        if len(t) == 0:
            continue
        
        t = t.replace('\t', ' ')
        tt = t.split(' ')
    
        while tt.count('') > 0:
            tt.remove('')

        if len(tt) == 0:
            continue

        if tt[1] == 'Unknown' or tt[1].count("Vent") == 1 \
           or tt[1].count("White-Matter") == 1 \
           or tt[1] == "Brain-Stem" \
           or tt[1].count("Accumbens") == 1:
            continue
    
        if tt[1].count("Left-") == 1:
            b_rois.append(tt[1].replace('Left-', 'lh_'))
        elif tt[1].count("Right-") == 1:
            b_rois.append(tt[1].replace('Right-', 'rh_'))

        b_ids.append(int(tt[0]))

    return (s_rois, s_ids, b_rois, b_ids)

if __name__ == '__main__':
    '''
    if len(sys.argv) < 3:
        print('Usage: gen_bips_aparc12_time_series.py sID imgMode [opts]')
        print('    imgMode = {fullspectrum | z_no_outliers_bandpassed}')
        print('    opts    = {-altaparc}')
        print('        -altaparc: specify an alternative aparc file, other than DATA/aparc12.nii.gz')
        sys.exit(0)
    '''

    parser = argparse.ArgumentParser(description= "Generate resting-state fMRI ROI time series based on the aparc12 parcellation")
    parser.add_argument("sID", type=str, help="Subject ID")
    parser.add_argument("imgMode", type=str, help="Image mode: {fullspectrum, z_no_outliers_bandpassed, z_no_outliers_bandpassed2, bpnrm, bpnrm2, bp2}")
    parser.add_argument("--altaparc", dest="altaparc", default="", \
                        help="Alternative aparc12 file name")
    parser.add_argument("--cuthead", dest="cuthead", default="", \
                        help="Remove the first specified number of  runs")
    parser.add_argument("--rebinarize", dest="rebinarize", action="store_true",\
                        help="Force re-generation of ROI masks with mri_binarize")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    
    sID = args.sID
    imgMode = args.imgMode
    altaparc = args.altaparc
    cuthead = args.cuthead
    rebinarize = args.rebinarize

    if len(cuthead) == 0:
        cuthead = 0
    else:
        cuthead = int(cuthead)
    if cuthead < 0:
        raise ValueError, "cuthead must be a positive integer."

    #sID = sys.argv[1]
    #imgMode = sys.argv[2]

    # ctabfn = os.path.join(FSDATA_dir, sID, 'label', 'aparc.annot.ctab')
    ctabfn = APARC12_TABLE

    '''
    [ids, rois] = read_ctab(ctabfn)
    sys.exit(0)
    # Expand the cortical rois into both hemispheres
    c_rois = []
    c_ids = []
    for (i0, hemi) in enumerate(HEMIS):
        for (i1, roi) in enumerate(rois):
            c_rois.append('%s_%s'%(hemi, roi))
            c_ids.append(1000 + 1000 * i0 + ids[i1])
    '''
   
    #hemi = sys.argv[2]

    print('sID = %s\n'%sID)
    #print('hemi = %s'%hemi)
    print("altaparc = %s\n"%altaparc)

    if len(altaparc) == 0:
        aparc_fn = os.path.join(DATA_dir, sID, 'aparc12.nii.gz')
    else:
        aparc_fn = altaparc
        print("INFO: Using altaparc = %s"%aparc_fn)

    if not os.path.isfile(aparc_fn):
        raise IOError, 'aparc file not found: %s'%aparc_fn
    print('aparc_fn = %s'%aparc_fn)

    if imgMode == 'z_no_outliers_bandpassed':
        resting_4d_fn = os.path.join(bips_resting_dir, sID, \
                                     'preproc', 'output', 'zscored', 'fwhm_5.0', \
                                     '%s_r00_z_no_outliers_bandpassed.nii.gz'%sID)
    elif imgMode == "z_no_outliers_bandpassed2":
        bips_resting_dir = bips_resting_dir_2
        resting_4d_fn = os.path.join(bips_resting_dir, sID, \
                                     'preproc', 'output', 'zscored', 'fwhm_0.0', \
                                     '%s_r00_z_no_outliers_bandpassed.nii.gz'%sID)
    elif imgMode == 'fullspectrum':
        resting_4d_fn = os.path.join(bips_resting_dir, sID, \
                                     'preproc', 'output', 'fullspectrum', 'fwhm_5.0', \
                                     '%s_r00_fullspectrum.nii'%sID)
    elif imgMode == "bpnrm":
        resting_4d_fn = os.path.join(bips_resting_dir, sID, \
                                     'preproc', 'output', 'bandpassed', \
                                     'fwhm_0.0', "%s_r00_bandpassed.nii"%(sID))
    elif imgMode == "bpnrm2":
        bips_resting_dir = bips_resting_dir_2
        resting_4d_fn = os.path.join(bips_resting_dir, sID, \
                                     "preproc", "output", "bandpassed", \
                                     "fwhm_0.0", "%s_r00_bandpassed.nii.gz"%(sID))
    elif imgMode == "bp2":
        bips_resting_dir = bips_resting_dir_2
        resting_4d_fn = os.path.join(bips_resting_dir, sID, \
                                     "preproc", "output", "bandpassed", \
                                     "fwhm_0.0", "%s_r00_bandpassed.nii.gz"%(sID))
    else:
        raise ValueError, 'Invalid imgMode: %s'%imgMode
        
    if not os.path.isfile(resting_4d_fn):
        raise IOError, '4D resting func file not found: %s'%resting_4d_fn
    print('resting_4d_fn = %s'%resting_4d_fn)

    if imgMode.endswith('2'):
        resting_mean_fn = os.path.join(bips_resting_dir, sID, \
                                       'preproc', 'mean', '%s_mean.nii.gz'%sID)
    else:
        resting_mean_fn = os.path.join(bips_resting_dir, sID, \
                                       'preproc', 'mean', '%s_mean.nii'%sID)
    if not os.path.isfile(resting_mean_fn):
        raise IOError, 'mean resting func file not found: %s'%resting_mean_fn
    print('resting_mean_fn = %s'%resting_mean_fn)
                                 
    bbreg_fsl_fn = os.path.join(bips_resting_dir, sID, \
                               'preproc', 'bbreg', '%s_register.mat'%sID)
    if not os.path.isfile(bbreg_fsl_fn):
        raise IOError, 'FSL-format bbreg file not found: %s'%bbreg_fsl_fn
    print('bbreg_fsl_fn = %s'%bbreg_fsl_fn)

    bbreg_inv_fsl_fn = os.path.join(bips_resting_dir, sID, \
                                    'preproc', 'bbreg', '%s_register_struct2func.mat'%sID)
    inv_xfm_cmd = 'convert_xfm -omat %s -inverse %s'%(bbreg_inv_fsl_fn, bbreg_fsl_fn)
    saydo(inv_xfm_cmd)
    
    #sys.exit(0)

    # Transform the aparc file to the resting-func space
    tmp_dir = os.path.join(ANALYSIS_DIR, sID, 'masks12')
    if not os.path.isdir(tmp_dir):
        os.system('mkdir -p %s'%tmp_dir)
        print('Created directory %s'%tmp_dir)
    else:
        #os.system('rm -r %s/*'%tmp_dir)
        #print('Cleaned directory %s'%tmp_dir)
        print("Directory already exists: %s"%tmp_dir)

    aparc_func_fn = os.path.join(tmp_dir, 'aparc_func.nii.gz')
    xfm_cmd = 'flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour'\
              %(aparc_fn, resting_mean_fn, bbreg_inv_fsl_fn, aparc_func_fn)
    saydo(xfm_cmd)

    # Generate the list of cortical and subcortical ROIs
    (rois, ids, sc_rois, sc_ids) = get_roi_ids(APARC12_TABLE, ASAP_SC_TABLE)
    
    c_rois = []
    c_ids = []
    for (i0, hemi) in enumerate(HEMIS):
        for (i1, roi) in enumerate(rois):
            c_rois.append('%s_%s'%(hemi, roi))
            c_ids.append(1000 + 1000 * i0 + int(ids[i1]))

    # 
    b_rois = c_rois + sc_rois
    b_ids = c_ids + sc_ids
    nc = len(c_rois) # Number of cortical ROIs
    
    # Determine the number of frames in the 4D resting fMRI file
    (stdout, stderr) = Popen(['mri_info', resting_4d_fn, '-P', '100'], \
                             stdout=PIPE).communicate()
    stdout = stdout.split('\n')
    bFound = False
    for t_line in stdout:
        if t_line.count('dimensions:') == 1:
            bFound = True
            nFrames = int(t_line.split(' ')[-1])

    if not bFound:
        raise ValueError, 'Unable to get fMRI series number of frames.'

    # Determine the outliers, if bpnrm mode is used
    if imgMode == "bpnrm":
        art_fn = os.path.join(bips_resting_dir, sID, 'preproc', 'art', \
                             'art._restingunwarped_outliers.txt')
    else:
        art_fn = os.path.join(bips_resting_dir, sID, 'preproc', 'art', \
                              'art._restingunwarped.nii_outliers.txt')

    if not os.path.isfile(art_fn):
        raise IOError, "Cannot find art outliers file: %s"%art_fn

    print("art_fn = %s"%art_fn)
    art_f = open(art_fn, 'r')
    art_txt = art_f.read().split('\n')
    art_f.close()
    outliers = []
    for t_line in art_txt:
        if len(t_line) > 0:
            outliers.append(int(t_line))

    if imgMode == "bpnrm" or imgMode == "bpnrm2":
        print("%d outliers found."%(len(outliers)))
        
        nFrames = nFrames - len(outliers)
        if len(outliers) > 0:
            print("nFrames --> %s"%(nFrames))
    
    # Process cut-head frames
    if cuthead > 0:
        if imgMode == "bpnrm" or imgMode == "bpnrm2":
            raise Exception, "Current, cuthead mode is not supported under bpnrm or bpnrm2 mode"

        chframes = []
        if cuthead >= nFrames:
            raise ValueError, "cuthead = %d >= nFrames = %d"%(cuthead, nFrames)
        for i0 in range(cuthead):
            chframes.append(i0)

        # Remove outlier time points that have already been removed by bips
        for olr in outliers:
            if chframes.count(olr) == 1:
                chframes.remove(olr)
                print("Removing outlier %d from chframes"%(olr))
                
        chframes0 = chframes
        chframes = []
        for i0 in range(len(chframes0)):
            chframes.append(i0)

        nFrames = nFrames - len(chframes)
        print("cuthead = %d: nFreames: %d --> %d"%(cuthead, nFrames + len(chframes), nFrames))
    else:
        chframes = []

    # Calculate the frame-by-frame in-brain mean intensity, for normalization
    if imgMode == "bpnrm" or imgMode == "bpnrm2":
        brainmean = np.zeros([nFrames])
        brainmask = os.path.join(bips_resting_dir, sID, 'preproc', 'mask', \
                                 "%s_brainmask.nii"%(sID))
        if not os.path.isfile(brainmask):
            raise IOError, "Canont find brain mask: %s"%(brainmask)
        
        masked_mean_cmd = "fslstats -t %s -k %s -m"%(resting_4d_fn, brainmask)
        (stdout, stderr) = Popen(masked_mean_cmd.split(' '), \
                                 stdout=PIPE, stderr=PIPE).communicate()
        meantxt = stdout.split('\n')
        cnt = 0
        for j0 in range(nFrames + len(outliers)):
            if outliers.count(j0) == 1:
                print("Skipping frame j0 = %d"%(j0))
                continue
            else:
                brainmean[cnt] = float(meantxt[j0])
                cnt = cnt + 1

        if len(np.nonzero(brainmean == 0)[0]) > 0:
            raise Exception, "Failed to calculate brain-wise intensity mean for all frames."

    nROIs = len(b_rois)
    bold_tab = np.array([[np.nan] * nFrames] * nROIs)
    
    for (i1, t_roi) in enumerate(b_rois):
        t_mask_fn = os.path.join(tmp_dir, 'mask_%s_rf.nii.gz'%t_roi)
        
        t_id = b_ids[i1]
        if os.path.isfile(t_mask_fn) and (not rebinarize):
            print('INFO: mask file already exists: %s'%t_mask_fn)
        else:
            binarize_cmd = 'mri_binarize --i %s --min %d --max %d --o %s'\
                           %(aparc_func_fn, t_id, t_id, t_mask_fn)
            saydo(binarize_cmd)

        #if t_id < 100:
        #    sys.exit(0)

        tmp_4d_fn = os.path.join(tmp_dir, 'tmp_4d_%s.nii.gz'%imgMode)
        multiply_cmd = 'fslmaths %s -mul %s %s'\
                       %(resting_4d_fn, t_mask_fn, tmp_4d_fn)
        saydo(multiply_cmd)

        #mean_cmd = 'fslstats -t %s -M'
        print('Extrating ROI-mean time course from ROI %s (%d) (%d / %d = %f%%)... \n'
              %(t_roi, t_id, i1, nROIs, float(i1) / float(nROIs) * 1e2))
        (stdout, stderr) = Popen(['fslstats', '-t', tmp_4d_fn, '-M'], 
                                 stdout=PIPE).communicate()
        t_bold_sig = stdout.split('\n')

        
        if (imgMode == "bpnrm" or imgMode == "bpnrm2") and len(outliers) > 0:
            t_bold_sig_0 = t_bold_sig
            t_bold_sig = []
            for (j0, t_val) in enumerate(t_bold_sig_0):
                if outliers.count(j0) == 1:
                    continue
                else:
                    t_bold_sig.append(t_val)
        elif len(chframes) > 0:
            t_bold_sig_0 = t_bold_sig
            t_bold_sig = []
            for (j0, t_val) in enumerate(t_bold_sig_0):
                if chframes.count(j0) == 1:
                    print("Cut head: skipping frame %d"%(j0))
                    continue
                else:
                    t_bold_sig.append(t_val)
                
        if imgMode == "bpnrm" or imgMode == "bpnrm2" : 
            # Do intensity normalization
            t_sig = []
            for (j0, t_val) in enumerate(t_bold_sig):
                if len(t_val) > 0:
                    t_sig.append(float(t_val))
            t_sig = np.array(t_sig)
            t_sig = t_sig / brainmean - 1
            t_sig = t_sig - np.mean(t_sig)

            for j1 in range(nFrames):
                bold_tab[i1][j1] = t_sig[j1]
        else:
            for j1 in range(nFrames):
                bold_tab[i1][j1] = float(t_bold_sig[j1])

    # Write immediate result to disk
    bold_tab_tmp_fn = os.path.join(tmp_dir, 'bold_tab_aparc12_tmp.pkl')
    fout = open(bold_tab_tmp_fn, 'wb')
    pickle.dump(bold_tab, fout)
    fout.close()
    print('bold_tab saved to (pickle) %s\n'%bold_tab_tmp_fn)

    # Compute the pairwise correlations 
    corr_tab = np.array([[np.nan] * nROIs] * nROIs)
    for i0 in range(nROIs):
        for i1 in range(nROIs):
            if (i1 <= i0):
                continue
            t_x = bold_tab[i0]
            t_y = bold_tab[i1]
            cc = np.corrcoef(t_x, t_y)
            corr_tab[i0][i1] = cc[0][1]

    #sys.exit(0)
            
    # Save final results to disk
    bips_resting_roi_corr = {'b_rois': b_rois, \
                             'b_ids': b_ids, \
                             'bold_tab': bold_tab, \
                             'corr_tab': corr_tab}
    out_fn = os.path.join(ANALYSIS_DIR, sID, 'roi_corr_aparc12.csc.%s.pkl'%imgMode)
    if cuthead > 0:
        out_fn = out_fn.replace(".pkl", ".ch%d.pkl"%(cuthead))
    fout = open(out_fn, 'wb')
    pickle.dump(bips_resting_roi_corr, fout)
    fout.close()
    print('Final results saved to (pickle) %s\n'%out_fn)
        
        
