#!/usr/bin/python

import os
import sys
import glob
import argparse
import numpy as np
from scipy.io import savemat
from subprocess import Popen, PIPE

from scai_utils import *
from aparc12 import get_aparc12_cort_rois

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
    ap = argparse.ArgumentParser(description="Generate single-hemisphere cortical WM connectivity matrix based on the aparc12_probtrackx results")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("hemi", help="hemisphere {lh, rh}")
    ap.add_argument("--speech", dest="bSpeech", action="store_true", \
                    help="Use the speech (sub-)network")
    ap.add_argument("--speechMode", type=str, default="", \
                    help="Speech network type (e.g., speech_PFS_lh)")
    ap.add_argument("--pt2", dest="bpt2", action="store_true", \
                    help="Use probtrackx2 results")

    if len(sys.argv) <= 1:
        ap.print_help()
        sys.exit(0)

    # Parse input arguments
    args = ap.parse_args()
    sID = args.sID
    hemi = args.hemi
    bpt2 = args.bpt2
    bSpeech = args.bSpeech
    speechMode = args.speechMode

    if bSpeech and len(speechMode) > 0:
        raise Exception, "Options --speech and --speechMode cannot be used together"
    if bSpeech:
        speechMode = "speech"
        
    if len(speechMode) > 0.05:
        if speechMode.count("_") !=3:
            raise Exception, \
                "Cannot find exactly 3 underlines in speechMode: %s" \
                % speechMode
        
    sm0 = speechMode.split("_")
    speechMode_noThr = "%s_%s_%s" % (sm0[0], sm0[1], sm0[2])

    if VALID_SPEECH_MODES.count(speechMode_noThr) == 0:
        raise Exception, "Unrecognized speechMode: %s" % speechMode

    if bpt2:
        TRACTS_RES_DIR = TRACTS_RES_DIR_PT2

    # Check sanity of input arguments
    if not (hemi == "lh" or hemi == "rh"):
        raise Exception, "Unrecognized hemisphere: %s" % hemi

    if len(speechMode) > 5:
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
    for t_roi in t_rois:
        h_rois.append(hemi + "_" + t_roi)

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

    # Rows: seeds; columns: targets
    for (i0, seedROI) in enumerate(h_rois):
        roiResDir = os.path.join(sResDir, seedROI)
        if bpt2:
            roiResDir += "_pt2"
        check_dir(roiResDir)
        
        fdtp = os.path.join(roiResDir, "fdt_paths.nii.gz")
        check_file(fdtp)

        fdtpn = os.path.join(roiResDir, "fdt_paths_norm.nii.gz")
        check_file(fdtpn)

        print("Processing seed ROI: %s" % seedROI)
        for (i1, targROI) in enumerate(h_rois):
            # Get mean
            mean_cmd = "fslstats %s -k %s -m" % (fdtp, diff_roi_masks[i1])
            (stdout, stderr) = Popen(mean_cmd.split(' '), \
                                     stdout=PIPE, stderr=PIPE)\
                               .communicate()
            if len(stderr) > 0:
                raise Exception, "Error occurred during %s" % mean_cmd
            t_val = float(stdout.split(' ')[0])
            connmat_mean[i0, i1] = t_val

            # Get mean normalized by seed size
            mean_norm_cmd = "fslstats %s -k %s -m" % (fdtpn, diff_roi_masks[i1])
            (stdout, stderr) = Popen(mean_norm_cmd.split(' '), \
                                     stdout=PIPE, stderr=PIPE)\
                               .communicate()
            if len(stderr) > 0:
                raise Exception, "Error occurred during %s" % mean_norm_cmd
            t_val = float(stdout.split(' ')[0])
            connmat_mean_norm[i0, i1] = t_val

            # Get median
            med_cmd = "fslstats %s -k %s -p 50" % (fdtp, diff_roi_masks[i1])
            (stdout, stderr) = Popen(med_cmd.split(' '), \
                                     stdout=PIPE, stderr=PIPE)\
                               .communicate()
            if len(stderr) > 0:
                raise Exception, "Error occurred during %s" % med_cmd
            t_val = float(stdout.split(' ')[0])
            connmat_median[i0, i1] = t_val

            # Get normalized median
            med_norm_cmd = "fslstats %s -k %s -p 50" \
                           % (fdtpn, diff_roi_masks[i1])
            (stdout, stderr) = Popen(med_norm_cmd.split(' '), \
                                     stdout=PIPE, stderr=PIPE)\
                               .communicate()
            if len(stderr) > 0:
                raise Exception, "Error occurred during %s" % med_norm_cmd
            t_val = float(stdout.split(' ')[0])
            connmat_median_norm[i0, i1] = t_val

    # Save results to mat file
    resMatFN = os.path.join(sResDir, "connmats.%s.mat" % hemi)    
    if bpt2:
        resMatFN = resMatFN.replace("connmats.", "connmats.pt2.")    
    if speechMode != "":
        resMatFN = resMatFN.replace("connmats.", "connmats.%s." % speechMode)

    savemat(resMatFN, \
            {"h_rois": h_rois, \
             "connmat_mean": connmat_mean, \
             "connmat_mean_norm": connmat_mean_norm, \
             "connmat_median": connmat_median, \
             "connmat_median_norm": connmat_median_norm})
    check_file(resMatFN)

    print("Connectivity matrices saved to mat file: %s" % resMatFN)
