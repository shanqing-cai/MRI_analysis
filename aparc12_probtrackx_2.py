#!/usr/bin/python

import os
import sys
import argparse
import glob
from scai_utils import *
from subprocess import Popen, PIPE

TRACULA_BASE = "/users/cais/STUT/analysis/dti2/tracula"
DATA_DIR = "/users/cais/STUT/DATA"
FSDATA_DIR = "/users/cais/STUT/FSDATA"
CTAB_FN = '/users/cais/STUT/slFRS17.ctab'
TRACTS_RES_DIR = "/users/cais/STUT/analysis/aparc12_tracts_2"
TRACTS_RES_DIR_PT2 = "/users/cais/STUT/analysis/aparc12_tracts_pt2"

TRACTSEG_BASE = "/users/cais/STUT/analysis/tractseg_aparc12"

PROJ_VOL_FWHM = 0
PROJ_SURF_FWHM = 0
TARG_SURF_FWHM = 5

SUBCORT_ROI_NAMES = ["Thalamus-Proper", "Caudate", "Putamen", "Pallidum"]

def get_roi_num(ctab_fn, roiHemi, roiName):
    # Read color table
    if not os.path.isfile(ctab_fn):
        raise Exception, "Cannot open color table (ctab): %s"%(ctab_fn)
    (ct_nums, ct_names) = read_ctab(ctab_fn)

    if not ct_names.count(roiName) == 1:
        raise Exception, "Not exactly one entry for ROI %s is found in%s"\
                         %(roiName, ctab_fn)

    roi_num = ct_nums[ct_names.index(roiName)]

    if roiHemi == "lh":
        roi_num += 1000
    elif roiHemi == "rh":
        roi_num += 2000
    else:
        raise Exception, "Unrecognized hemisphere name: %s"%roiHemi

    return roi_num


if __name__ == "__main__":
    ap = argparse.ArgumentParser("Run probtrackx with various options")
    ap.add_argument("sID", type=str, help="Subject ID")
    ap.add_argument("aparc12_roiName", type=str, help="Name of ROI in the aparc12 paradigm")
    ap.add_argument("--stoplabel", type=str, dest="stoplabel", default="", \
                    help="Stop mask label on the surface. Format: subjID,labelFN. Example: fsaverage,/users/cais/STUT/analysis/aparc12_tracts/labels/lh_IFS_rsfc_lh_vIFo_FWHM5.label")
    ap.add_argument("--ccavoid", dest="bCCAvoid", action="store_true", \
                    help="Use corpus callosum stopmask")
    ap.add_argument("--ccwaypoint", dest="bCCWaypoint", action="store_true", 
                    help="Use corpus callosum waypoint (incompatible with --ccavoid; mainly for cross-hemisphere tracking")
    ap.add_argument("--ww", dest="bWMWaypoint", action="store_true", \
                    help="Use ipsilateral white-matter as waypoint")
    ap.add_argument("--noProbtrackx", dest="bNoProbtrackx", \
                    action="store_true", \
                    help="Don't do the actual time-consuming tracking.")
    ap.add_argument("--pt2", dest="bpt2", action="store_true", \
                    help="Use probtrackx2 [FALSE]")
    ap.add_argument("--noCleanUp", dest="bNoCleanUp", action="store_true", \
                    help="Don't perform clean up.")
    ap.add_argument("--outdir", type=str, dest="uOutDir", default="", \
                    help="User specified output directory")
    ap.add_argument("--noSurfProj", dest="bNoSurfProj", action="store_true", \
                    help="Do not perform projection to surface")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
    
    args = ap.parse_args()
    sID = args.sID
    aparc12_roiName = args.aparc12_roiName
    stoplabel = args.stoplabel
    bCCAvoid = args.bCCAvoid
    bCCWaypoint = args.bCCWaypoint
    bWMWaypoint = args.bWMWaypoint
    bNoProbtrackx = args.bNoProbtrackx
    bNoCleanUp = args.bNoCleanUp
    bpt2 = args.bpt2
    uOutDir = args.uOutDir
    bNoSurfProj = args.bNoSurfProj

    """
    if bCCAvoid and stoplabel == "":
        raise Exception, "--ccavoid must be used together with --stoplabel"
    """

    # Locate the subject data directory
    sDataDir = os.path.join(DATA_DIR, sID)
    check_dir(sDataDir)
    
    aparc12img = os.path.join(sDataDir, "aparc12.nii.gz")
    check_file(aparc12img)

    if bpt2:
        TRACTS_RES_DIR = TRACTS_RES_DIR_PT2

    # Get the ROI mask in the diffision space
    sResDir = os.path.join(TRACTS_RES_DIR, sID)
    check_dir(sResDir, bCreate=True)

    # Prepare CC (corpus callosum) mask (optional)
    if bCCAvoid:
        ccMask = os.path.join(FSDATA_DIR, sID, "mri", "ccMask_diff.nii.gz")
        check_file(ccMask)

    # Parse stoplabel (optional)
    if stoplabel != "":
        stopLblVolFN = os.path.join(sResDir, \
                                    "aparc12_%s.diff.nii.gz" % stoplabel)
        check_file(stopLblVolFN)
        """
        stoplbl_src_id = stoplabel.split(',')[0]
        stoplbl_src_fn = stoplabel.split(',')[1]
        check_file(stoplbl_src_fn)

        (foo, stoplbl_basefn) = os.path.split(stoplbl_src_fn)
        if stoplbl_basefn.startswith("lh"):
            stopLblHemi = "lh"
        elif stoplbl_basefn.startswith("rh"):
            stopLblHemi = "rh"
        else:
            raise Exception, "Cannot determine the hemisphere of input stop label: %s"%stoplbl_src_fn
        
        # Convert the label on the source surface to a stop mask in the destination volume
        sStopLblDir = os.path.join(sResDir, "stoplabels")
        check_dir(sStopLblDir, bCreate=True)
        sLblFN = os.path.join(sStopLblDir, stoplbl_basefn)
        mapcmd = "mri_label2label --srclabel %s --srcsubject %s "\
                 %(stoplbl_src_fn, stoplbl_src_id) + \
                 "--trglabel %s --trgsubject %s "%(sLblFN, sID) + \
                 "--regmethod surface --hemi %s"%(stopLblHemi)        
        os.system(mapcmd)
        check_file(sLblFN)

        stopLblVolFN = os.path.join(sStopLblDir, \
                                    stoplbl_basefn.replace(".label", ".nii.gz"))
        proj_cmd = "mri_label2vol --subject %s --hemi %s --label %s "\
                   %(sID, stopLblHemi, sLblFN) + \
                   "--proj frac 0 1 0.1 --identity --temp %s --o %s"\
                   %(aparc12img, stopLblVolFN)
        os.system(proj_cmd)
        check_file(stopLblVolFN)
        """
        
    # Determine whether the ROI is a subcortical ROI
    if not (aparc12_roiName.startswith("lh_") or \
            aparc12_roiName.startswith("rh_")):
        raise Exception, "roiName must start with lh_ or rh_"

    bSC = False
    for (i0, scn) in enumerate(SUBCORT_ROI_NAMES):
        if aparc12_roiName.count(scn) == 1:
            bSC = True
            break

    if bSC:
        print("INFO: detected subcortical ROI name")

        # TODO: segB and so on
        if not aparc12_roiName.split(".")[-1].startswith("segA"):
            raise Exception, "Subcoritcal ROI name: cannot find seg number"
        
        scSegNumStr = aparc12_roiName.split(".")[-1].split("segA")[-1]        

        if scSegNumStr.count('-') == 0:
            scSegNum = [-1]
            scSegNum[0] = int(scSegNumStr)

            print("INFO: subcortical segmentation name = %d" % scSegNum[0])
        elif scSegNumStr.count('-') == 1:
            scSegNum = [-1, -1]
            scSegNum[0] = int(scSegNumStr.split("-")[0])
            scSegNum[1] = int(scSegNumStr.split("-")[1])

            print("INFO: subcortical segmentation name = %d - %d" % \
                  (scSegNum[0], scSegNum[1]))
        else:
            raise Exception, "Unrecognized scSegNumStr: %s" % scSegNumStr
        

    # === Prepare WM waypoint mask (optional) ===
    sTracDir = os.path.join(TRACULA_BASE, sID)
    check_dir(sTracDir)

    FA_img = os.path.join(sTracDir, "dmri", "dtifit_FA.nii.gz")
    check_file(FA_img)
    
    if bWMWaypoint:
        if bSC:
            raise Exception, "White-matter waypoint (WM) is not currently supported under bSC = True"
        hemi = aparc12_roiName.split('_')[0]
        ipsiWMMask = os.path.join(FSDATA_DIR, sID, \
                                  "mri", "wmMask_%s_diff.nii.gz" % hemi)
        if os.path.isfile(ipsiWMMask):
            print("INFO: Ipsilateral white-matter waypoint mask already exisits: %s" % ipsiWMMask)
        else:
            # Generate it 
            ipsiWMStruct1 = os.path.join(FSDATA_DIR, sID, \
                                        "mri", "wmMask_%s_1.nii.gz" % hemi)
            ipsiWMStruct2 = os.path.join(FSDATA_DIR, sID, \
                                        "mri", "wmMask_%s_2.nii.gz" % hemi)
            ipsiWMStruct = os.path.join(FSDATA_DIR, sID, \
                                        "mri", "wmMask_%s.nii.gz" % hemi)

            if hemi == "lh":
                idxRange = [3000, 3999]
                idxExtra = 5001
            else:
                idxRange = [4000, 4999]
                idxExtra = 5002
                
            bin_cmd = "mri_binarize --i %s --min %d --max %d --o %s" % \
                      (aparc12img, idxRange[0], idxRange[1], ipsiWMStruct1)
            os.system("rm -f %s" % ipsiWMStruct1)
            saydo(bin_cmd)
            check_file(ipsiWMStruct1)
            
            bin_cmd = "mri_binarize --i %s --min %d --max %d --o %s" % \
                      (aparc12img, idxExtra, idxExtra, ipsiWMStruct2)
            os.system("rm -f %s" % ipsiWMStruct2)
            saydo(bin_cmd)
            check_file(ipsiWMStruct2)

            add_cmd = "fslmaths %s -add %s %s" % \
                      (ipsiWMStruct1, ipsiWMStruct2, ipsiWMStruct)
            os.system("rm -f %s" % ipsiWMStruct)
            saydo(add_cmd)
            check_file(ipsiWMStruct)

            # == Transform to diffusion space == #
            a2d = os.path.join(sTracDir, "dmri", "xfms", "a2d.mat")
            check_file(a2d)

            xfm_cmd = "flirt -in %s -ref %s -applyxfm -init %s -out %s" % \
                      (ipsiWMStruct, FA_img, a2d, ipsiWMMask) + \
                      " -interp nearestneighbour"
            os.system("rm -f %s" % ipsiWMMask)
            saydo(xfm_cmd)
            check_file(ipsiWMMask)
                

    # === Create the seed mask volume === %    
    mask_img = os.path.join(sTracDir, "dmri", "nodif_brain_mask.nii.gz")
    check_file(mask_img)

    if not bSC:
        # Figure out the ROI code for the ROI
        roiHemi = aparc12_roiName.split('_')[0]
        roiName = aparc12_roiName.split('_')[1]
    
        roi_num = get_roi_num(CTAB_FN, roiHemi, roiName)

        anat2diff_xfm_mat = os.path.join(sTracDir, "dmri", "xfms", \
                                     "anatorig2diff.bbr.mat")
#    anat2diff_xfm_mat = os.path.join(sTracDir, "dmri", "xfms", \
#                                     "diff2anatorig.bbr.mat")
        check_file(anat2diff_xfm_mat)

        bin_out = os.path.join(sResDir, "aparc12_%s.nii.gz"%aparc12_roiName)
        xfm_out = os.path.join(sResDir, \
                               "aparc12_%s.diff.nii.gz"%aparc12_roiName)

        if not os.path.isfile(xfm_out):
            bin_cmd = "mri_binarize --i %s --match %d --o %s" \
                      % (aparc12img, roi_num, bin_out)
            saydo(bin_cmd)
            check_file(bin_out)

            xfm_cmd = "flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour" \
                %(bin_out, FA_img, anat2diff_xfm_mat, xfm_out)
            saydo(xfm_cmd)
            check_file(xfm_out)

        """
    if stoplabel != "":
        xfm_stopLbl_out = stopLblVolFN.replace(".nii.gz", ".diff.nii.gz")

#        if bCCAvoid:
#            xfm_stopLbl_out = xfm_stopLbl_out.replace(".diff.nii.gz", \
#                                                      ".ccAvoid.diff.nii.gz")

        xfm_stopLbl_cmd = "flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour"\
            %(stopLblVolFN, FA_img, anat2diff_xfm_mat, xfm_stopLbl_out)
        os.system(xfm_stopLbl_cmd)
        check_file(xfm_stopLbl_out)

        #if bCCStop:
        #    add_cmd = "fslmaths %s -add %s %s" % \
        #              (xfm_stopLbl_out, ccMask, xfm_stopLbl_out)
        #    saydo(add_cmd)
        """
        
        seedMaskFN = xfm_out
    else:
        scROIName = aparc12_roiName.split(".segA")[0].replace("h_", "h.")
        tractSegFN = os.path.join(TRACTSEG_BASE, sID, scROIName + "_ccStop", \
                                  "tract_regionParc.nii.gz")
        check_file(tractSegFN)

        if len(scSegNum) == 1:
            tractSegSeedMask = tractSegFN.replace(".nii.gz", \
                                              "_seg%d.nii.gz" % scSegNum[0])

            bin_cmd = "mri_binarize --i %s --match %d --o %s" \
                      % (tractSegFN, scSegNum[0], tractSegSeedMask)
        else:
            tractSegSeedMask = tractSegFN.replace(".nii.gz", \
                                              "_seg%d-%d.nii.gz" % \
                                              (scSegNum[0], scSegNum[1]))

            bin_cmd = "mri_binarize --i %s --min %d --max %d --o %s" \
                   % (tractSegFN, scSegNum[0], scSegNum[1], tractSegSeedMask)
        
        delete_file_if_exists(tractSegSeedMask)
        saydo(bin_cmd)
        check_file(tractSegSeedMask)

        seedMaskFN = tractSegSeedMask
        
    # Locate the bedpostx base dir
    bpbase = os.path.join(TRACULA_BASE, sID, "dmri.bedpostX")
    check_dir(bpbase)
    bpbase = os.path.join(bpbase, "merged")

    mfns = glob.glob(bpbase + "*")
    if len(mfns) < 6:
        raise Exception, "It appears that the bedpostx data are incomplete at: %s"%(bpbase + "*")

    # Provide command line for checking the coregistration of DTI and anatomical volumes:
    check_reg_cmd = "tkregister2 --mov %s --targ %s --identity --reg identity.dat" \
                    % (seedMaskFN, FA_img)
    print("Command for checking the diffusion-to-anatomical coregsitration:\n\t%s\n" % (check_reg_cmd))

    # Run probtrackx
    if stoplabel == "":
        outdir = os.path.join(sResDir, aparc12_roiName)
    else:
        """
        outdir = os.path.join(sResDir, \
                              aparc12_roiName + "-" + \
                              stoplbl_basefn.replace(".label", ""))
        """
        outdir = os.path.join(sResDir, \
                              aparc12_roiName + "-" + \
                              stoplabel)

    if bCCAvoid:
        if not bWMWaypoint:
            outdir += "_ccAvoid"
        else:
            outdir += "_caww"

    if bpt2:
        outdir += "_pt2"

    if len(uOutDir) > 0:
        outdir = uOutDir

    if not bpt2:
        cmd = 'probtrackx --mode=seedmask -x ' + seedMaskFN + \
              ' -s ' + bpbase + \
              ' -m ' + mask_img + ' -l -c 0.2 -S 2000 --steplength=0.5' + \
              ' -P 5000 ' + \
              ' --forcedir --opd --pd --dir=' + outdir
    else:
        cmd = 'probtrackx2 -x ' + seedMaskFN + \
              ' -s ' + bpbase + \
              ' -m ' + mask_img + ' -l -c 0.2 -S 2000 --steplength=0.5' + \
              ' -P 5000 ' + \
              ' --forcedir --opd --pd --dir=' + outdir

    # === Destination label === #
    if stoplabel != "":
        """
        wptext = xfm_stopLbl_out + '\n'        
        wpfn = stopLblVolFN.replace(".nii.gz", ".txt")
        wpf = open(wpfn, "w")
        wpf.write(wptext)
        wpf.close()
        check_file(wpfn)
        cmd += " --stop=%s --waypoints=%s"%(xfm_stopLbl_out, wpfn)
        """

        wptext = stopLblVolFN + "\n"
        wpfn = stopLblVolFN.replace(".nii.gz", ".txt")
        wpf = open(wpfn, "w")
        wpf.write(wptext)
        wpf.close()
        check_file(wpfn)
        cmd += " --stop=%s --waypoints=%s"%(stopLblVolFN, wpfn)
        #cmd += " --waypoints=%s --stop=%s"%(wpfn, stopLblVolFN)

        """
        check_targ_reg_cmd = "tkregister2 --mov %s --targ %s --identity --reg identity.dat" \
                             % (xfm_stopLbl_out, FA_img)
        """
        check_targ_reg_cmd = "tkregister2 --mov %s --targ %s --identity --reg identity.dat" \
                             % (stopLblVolFN, FA_img)

    # === Ipsilateral WM waypoint === #
    if bWMWaypoint:
        if stoplabel != "":
            raise Exception, "Current version does not support simultaneous bWMWaypoint and stoplabel"

        wptext = ipsiWMMask + "\n"
        wpfn = ipsiWMMask.replace(".nii.gz", ".txt")
        wpf = open(wpfn, "w")
        wpf.write(wptext)
        wpf.close()

        check_file(wpfn)
        cmd += " --waypoints=%s" % (wpfn)
        
    if bCCAvoid:
        cmd += " --avoid=%s" % ccMask

    if not bNoProbtrackx:
        saydo(cmd)

    check_dir(outdir)

    # Get the sizes of the seed (and the stop mask) and write the info to the results directory
    (sout, serr) = Popen(['fslstats', seedMaskFN, '-V'], 
                         stdout=PIPE, stderr=PIPE).communicate()
    seed_size_fn = os.path.join(outdir, 'seed_size.txt')
    seed_size_f = open(seed_size_fn, 'w')
    seed_size_f.write(sout)
    seed_size_f.close()
    check_file(seed_size_fn)

    seed_nvox = int(sout.split(' ')[0])

    # Normalize the fdt_path file
    fdt_paths_fn = os.path.join(outdir, "fdt_paths.nii.gz")
    check_file(fdt_paths_fn)
    
    fdt_paths_norm_fn = os.path.join(outdir, "fdt_paths_norm.nii.gz")
    norm_cmd = "fslmaths -dt float %s -div %d %s -odt float" % \
               (fdt_paths_fn, seed_nvox, fdt_paths_norm_fn)
    saydo(norm_cmd)
    check_file(fdt_paths_norm_fn)

    if stoplabel != "":
        """
        (sout, serr) = Popen(['fslstats', xfm_stopLbl_out, '-V'], 
                              stdout=PIPE, stderr=PIPE).communicate()
        """
        (sout, serr) = Popen(['fslstats', stopLblVolFN, '-V'], 
                              stdout=PIPE, stderr=PIPE).communicate()
        targ_size_fn = os.path.join(outdir, 'targ_size.txt')
        targ_size_f = open(targ_size_fn, 'w')
        targ_size_f.write(sout)
        targ_size_f.close()
        check_file(targ_size_fn)

    # Project to surface and then to fsaverage
    
#    proj_reg = os.path.join(sTracDir, "dmri", "xfms", \
#                            "anatorig2diff.bbr.dat")
    proj_reg = os.path.join(sTracDir, "dmri", "xfms", \
                            "diff2anatorig.bbr.dat")
    #check_file(proj_reg)

    if not os.path.isfile(proj_reg):
        # Convert the fsl-style mat to fs-style dat
        proj_reg_fsl = os.path.join(sTracDir, "dmri", "xfms", \
                                    "diff2anatorig.bbr.mat")
        check_file(proj_reg_fsl)
        tkr_cvt_cmd = "tkregister2 --mov %s --targ %s --fsl %s " % \
                      (FA_img, aparc12img, proj_reg_fsl) + \
                      "--s %s --reg %s --noedit" % (sID, proj_reg)

        saydo(tkr_cvt_cmd)
        check_file(proj_reg)
            
    hemis = ["lh", "rh"]
    fsav_surf_view_cmd = {}
    for hemi in hemis:
        proj_out = os.path.join(outdir, "fdt_paths_norm.%s.mgz" % hemi)
        proj_cmd = "mri_vol2surf --mov %s --reg %s --o %s --hemi %s --trgsubject %s --projfrac 0 --fwhm %d --surf-fwhm %d --noreshape --interp trilin" \
                   % (fdt_paths_norm_fn, proj_reg, proj_out, hemi, \
                      sID, PROJ_VOL_FWHM, PROJ_SURF_FWHM)

        if not bNoSurfProj:
            saydo(proj_cmd)
            check_file(proj_out)

        fsav_fn = os.path.join(outdir, "fdt_paths_norm.fsav.%s.mgz" % hemi)
        preproc_cmd = "mris_preproc --s %s --hemi %s --is %s --target fsaverage --fwhm %d  --out %s" \
                      % (sID, hemi, proj_out, TARG_SURF_FWHM, fsav_fn)

        if not bNoSurfProj:
            saydo(preproc_cmd)
            check_file(fsav_fn)

        fsav_surf_view_cmd[hemi] = "tksurfer fsaverage %s inflated -gray -overlay %s -fmid 100" % (hemi, fsav_fn)

    if not bNoCleanUp:
        if not bSC:
            saydo("rm -f %s"%bin_out)
            #saydo("rm -f %s"%xfm_out)

    # Print viewing command    
    check_file(fdt_paths_fn)
    viewCmd = "freeview %s %s:colormap=jet"%(FA_img, fdt_paths_fn)
    print("To view the result, do:\n\t%s\n"%viewCmd)

    for hemi in hemis:
        print("To view the results on fsaverage %s surface, do:\n\t%s\n" \
              % (hemi, fsav_surf_view_cmd[hemi]))

    print("\nCommand for checking the diffusion-to-anatomical coregsitration:\n\t%s" % (check_reg_cmd))

    if stoplabel != "":
        print("\nCommand for checking the diffusion-to-anatomical coregsitration for the target:\n\t%s" % (check_targ_reg_cmd))

    # Write the view commands to file: viewcmds.txt
    viewcmds_fn = os.path.join(outdir, "viewcmds.txt")
    viewcmds_f = open(viewcmds_fn, 'w')
    
    viewcmds_f.write(check_reg_cmd + '\n\n')
    if stoplabel != "":
        viewcmds_f.write(check_targ_reg_cmd + '\n\n')

    viewcmds_f.write(viewCmd + '\n\n')

    for hemi in hemis:
        viewcmds_f.write(fsav_surf_view_cmd[hemi] + "\n\n")

    viewcmds_f.close()
    check_file(viewcmds_fn)

    
