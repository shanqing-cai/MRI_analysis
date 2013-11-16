#!/usr/bin/python

import os
import sys
import glob
import argparse
import numpy as np
from scipy.io import savemat
from subprocess import Popen, PIPE
import nibabel as nb

from scai_utils import *

DATA_DIR = "/users/cais/STUT/DATA"
TRACULA_BASE = "/users/cais/STUT/analysis/dti2/tracula"
TRACTS_RES_DIR = "/users/cais/STUT/analysis/aparc12_tracts_2"
TRACTS_RES_DIR_PT2 = "/users/cais/STUT/analysis/aparc12_tracts_pt2"

CTAB_FN_CORTICAL = '/users/cais/STUT/slFRS17.ctab'
CTAB_FN_SUBCORTICAL = '/software/atlas/ASAP_subcortical_labels.txt'

VALID_SPEECH_MODES = ['', 'speech', 'speech_PFS_lh', 'speech_PFS_rh', \
                      'speech_PWS_lh', 'speech_PWS_rh', \
                      'speech_2g_lh', 'speech_2g_rh']

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate single-hemisphere (lh/rh) or cross-hemisphere (xh) cortical WM connectivity matrix based on the aparc12_probtrackx results")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("hemi", help="hemisphere {lh - left hemisphere, rh - right hemisphere, xh - cross-hemisphere}")
    ap.add_argument("--speech", dest="bSpeech", action="store_true", \
                    help="Use the speech (sub-)network")
    ap.add_argument("--speechMode", type=str, default="", \
                    help="Speech network type (e.g., speech_PFS_lh)")
    ap.add_argument("--caww", dest="bCAWW", action="store_true", \
                    help="Corpus-callosum avoidance and ipsilateral WM waypoint mask")
    ap.add_argument("--cw", dest="bCW", action="store_true", \
                    help="Corpos-callosum waypoint (for cross-hemisphere tracking (not compatible with option --caww)")
    ap.add_argument("--pt2", dest="bpt2", action="store_true", \
                    help="Use probtrackx2 results")
    ap.add_argument("--oldver", dest="bOldVer", action="store_true",
                    help="Use the old version of ROI names")

    if len(sys.argv) <= 1:
        ap.print_help()
        sys.exit(0)

    # Parse input arguments
    args = ap.parse_args()
    sID = args.sID
    hemi = args.hemi
    bpt2 = args.bpt2
    bCAWW = args.bCAWW
    bCW = args.bCW
    bSpeech = args.bSpeech
    speechMode = args.speechMode

    if bCAWW and bCW:
        raise Exeption, "Using incompatible options: --caww and --cw"

    if args.bOldVer:
        from aparc12_oldVer import get_aparc12_cort_rois
    else:
        from aparc12 import get_aparc12_cort_rois

    if bSpeech and len(speechMode) > 0:
        raise Exception, "Options --speech and --speechMode cannot be used together"
    if bSpeech:
        speechMode = "speech"

    if len(speechMode) > 6:
        if speechMode.count("_") !=3:
            raise Exception, \
                "Cannot find exactly 3 underlines in speechMode: %s" \
                % speechMode
        
        sm0 = speechMode.split("_")
        speechMode_noThr = "%s_%s_%s" % (sm0[0], sm0[1], sm0[2])
    else:
        speechMode_noThr = ""

    if VALID_SPEECH_MODES.count(speechMode_noThr) == 0:
        raise Exception, "Unrecognized speechMode: %s" % speechMode

    if bpt2:
        TRACTS_RES_DIR = TRACTS_RES_DIR_PT2
    check_dir(TRACTS_RES_DIR)

    # Check sanity of input arguments
    if not (hemi == "lh" or hemi == "rh" or hemi == "xh"):
        raise Exception, "Unrecognized hemisphere: %s" % hemi

    if len(speechMode) > 6:
        if hemi != speechMode_noThr[-2 :]:
            raise Exception, "Mismatch between hemi=%s and speechMode=%s" \
                             % (hemi, speechMode)

    # Read cortical ROI list
    if speechMode == "speech":
        t_rois = get_aparc12_cort_rois("all", bSpeech=True)
    else:
        t_rois = get_aparc12_cort_rois("all", bSpeech=speechMode)
    t_rois.sort()
    
    h_rois = []
    if hemi == "lh" or hemi == "rh":    
        for t_roi in t_rois:
            h_rois.append(hemi + "_" + t_roi)
    elif hemi == "xh":
        hemis = ["lh", "rh"]
        for t_hemi in hemis:
            for t_roi in t_rois:
                h_rois.append(t_hemi + "_" + t_roi)

    # Preparation: check directories
    sDataDir = os.path.join(DATA_DIR, sID)
    check_dir(sDataDir)

    sTracDir = os.path.join(TRACULA_BASE, sID)
    check_dir(sTracDir)

    sResDir = os.path.join(TRACTS_RES_DIR, sID)
    check_dir(sResDir)

    # Preparation: check the existence of all diffusion-space ROI masks
    diff_roi_masks = []
    for h_roi in h_rois:
        diff_roi_mask = os.path.join(sResDir, \
                                    "aparc12_%s.diff.nii.gz" % h_roi)
        check_file(diff_roi_mask)
        diff_roi_masks.append(diff_roi_mask)
    nROIs = len(h_rois)

    # Iterate through all seed ROIs and build up the connectivity matrix
    connmat_mean = np.zeros([nROIs, nROIs])
    connmat_median = np.zeros([nROIs, nROIs])
    connmat_mean_norm = np.zeros([nROIs, nROIs])
    connmat_median_norm = np.zeros([nROIs, nROIs])

    #=== Load ROI masks ===#
    info_log("Loading %d ROI masks" % len(diff_roi_masks))

    mask_shapes = []
    nzIdx = []
    for (i0, roi_mask) in enumerate(diff_roi_masks):
        check_file(roi_mask)
        
        t_img = nb.load(roi_mask)
        t_img_dat = t_img.get_data()

        mask_shapes.append(np.shape(t_img_dat))
            
        t_img_dat = np.ndarray.flatten(t_img_dat)
            
        nzIdx.append(np.nonzero(t_img_dat)[0])

    if len(np.unique(mask_shapes)) != 1:
        error_log("Non-unique matrix size among the mask files", logFN=logFN)
    imgShape = np.unique(mask_shapes)[0]

    #=== Get connectivity data ===#
    # Rows: seeds; columns: targets
    for (i0, seedROI) in enumerate(h_rois):
        roiResDir = os.path.join(sResDir, seedROI)
        if bCAWW:
            roiResDir += "_caww"
        elif bCW:
            roiResDir += "_cw"

        if bpt2:
            roiResDir += "_pt2"

        check_dir(roiResDir)
        
        fdtp = os.path.join(roiResDir, "fdt_paths.nii.gz")
        check_file(fdtp)

        fdtpn = os.path.join(roiResDir, "fdt_paths_norm.nii.gz")
        check_file(fdtpn)

        info_log("Processing seed ROI: %s" % seedROI)

        # Load fdtp image 
        t_img = nb.load(fdtp)
        t_img_dat = t_img.get_data()
        assert(list(np.shape(t_img_dat)) == list(imgShape))
        img_dat = np.ndarray.flatten(t_img_dat)

        # Load fdtpn image
        t_img = nb.load(fdtpn)
        t_img_dat = t_img.get_data()
        assert(list(np.shape(t_img_dat)) == list(imgShape))
        img_dat_n = np.ndarray.flatten(t_img_dat)
        
        for (i1, targROI) in enumerate(h_rois):
            if hemi == "xh" and seedROI[:3] == targROI[:3]:
                continue # xh: Omit same-hemisphere projections

            connmat_mean[i0, i1] = np.mean(img_dat[nzIdx[i1]])
            connmat_mean_norm[i0, i1] = np.mean(img_dat_n[nzIdx[i1]])
            connmat_median[i0, i1] = np.median(img_dat[nzIdx[i1]])
            connmat_median_norm[i0, i1] = np.median(img_dat_n[nzIdx[i1]])

    # Save results to mat file
    resMatFN = os.path.join(sResDir, "connmats.%s.mat" % hemi)
    if bpt2:
        resMatFN = resMatFN.replace("connmats.", "connmats.pt2.")
    if bCAWW:
        resMatFN = resMatFN.replace("connmats.", "connmats.caww.")
    if bCW:
        resMatFN = resMatFN.replace("connmats.", "connmats.cw.")
    if speechMode != "":
        resMatFN = resMatFN.replace("connmats.", "connmats.%s." % speechMode)

    savemat(resMatFN, \
            {"h_rois": h_rois, \
             "connmat_mean": connmat_mean, \
             "connmat_mean_norm": connmat_mean_norm, \
             "connmat_median": connmat_median, \
             "connmat_median_norm": connmat_median_norm})
    check_file(resMatFN)

    info_log("Connectivity matrices saved to mat file: %s" % resMatFN)
