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
TRACT_SEG_BASE = "/users/cais/STUT/analysis/tractseg_aparc12"
TRACTS_RES_DIR = "/users/cais/STUT/analysis/aparc12_tracts_2"
TRACTS_RES_DIR_PT = "/users/cais/STUT/analysis/aparc12_tracts_2"
TRACTS_RES_DIR_PT2 = "/users/cais/STUT/analysis/aparc12_tracts_pt2"

CTAB_FN_CORTICAL = '/users/cais/STUT/slFRS17.ctab'
CTAB_FN_SUBCORTICAL = '/software/atlas/ASAP_subcortical_labels.txt'

SEG_TYPE_NUMSEGS = {"A": 7, "A2_5": 1, "A2_4": 1}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate single-hemisphere WM connectivity vector from specified subcortical seed with tract_seg to aparc12 cortical ROIs")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("hemi", help="hemisphere {lh, rh}")
    ap.add_argument("scROI", help="Subcortical ROI name (e.g., Putamen)")
    ap.add_argument("segType", help="Tract-seg type (e.g., A, the default)")
    # ap.add_argument("segNum", type=int, help="Tract-seg number (e.g., 1)")
    ap.add_argument("--ccstop", dest="bCCStop", action="store_true", \
                    help="Use corpus callosum stop mask results");
    ap.add_argument("--speech", dest="bSpeech", action="store_true", \
                    help="Use the speech (sub-)network")
    ap.add_argument("--pt2", dest="bpt2", action="store_true", \
                    help="Use probtrackx2 results")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    # Parse input arguments
    args = ap.parse_args()
    
    sID = args.sID
    hemi = args.hemi
    scROI = args.scROI
    segType = args.segType
    #segNum = args.segNum
    
    bCCStop = args.bCCStop
    bSpeech = args.bSpeech
    bpt2 = args.bpt2

    if bpt2:
        TRACTS_RES_DIR = TRACTS_RES_DIR_PT2

    # Check sanity of input arguments
    if not (hemi == "lh" or hemi == "rh"):
        raise Exception, "Unrecognized hemisphere: %s" % hemi

    if segType.count("-") == 1:
        segType1 = segType.replace("-", "_")
        if SEG_TYPE_NUMSEGS.keys().count(segType1) != 1:
            raise Exception, "Unrecognized segType: %s" % segType1
        nSegs = SEG_TYPE_NUMSEGS[segType1]
    else:
        if SEG_TYPE_NUMSEGS.keys().count(segType) != 1:
            raise Exception, "Unrecognized segType: %s" % segType
        nSegs = SEG_TYPE_NUMSEGS[segType]

    # Read cortical ROI list
    t_rois = list(get_aparc12_cort_rois(bSpeech=bSpeech))
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
    
    sTractSegDir = os.path.join(TRACT_SEG_BASE, sID)
    check_dir(sTractSegDir)

    # Preparation: check the existence of all diffusion-space ROI masks
    diff_roi_masks = []
    for h_roi in h_rois:
        diff_roi_mask = os.path.join(TRACTS_RES_DIR_PT, sID, \
                                    "aparc12_%s.diff.nii.gz" % h_roi)
        check_file(diff_roi_mask)
        diff_roi_masks.append(diff_roi_mask)
    nROIs = len(h_rois)

    # Iterate through all seed ROIs and build up the connectivity matrix
    connmat_mean = np.zeros([nSegs, nROIs])
    connmat_median = np.zeros([nSegs, nROIs])
    connmat_mean_norm = np.zeros([nSegs, nROIs])
    connmat_median_norm = np.zeros([nSegs, nROIs])

    # 
    if not bCCStop:
        segResDir = os.path.join(sTractSegDir, hemi + "." + scROI)
    else:
        segResDir = os.path.join(sTractSegDir, hemi + "." + scROI + "_ccStop")

    check_dir(segResDir)
    
    # Rows: seeds; columns: targets
    for i0 in range(nSegs):
        segn = i0 + 1

        if nSegs == 1:
            roiResDir = os.path.join(sResDir, \
                       "%s_%s.seg%s" % (hemi, scROI, segType))
        else:
            roiResDir = os.path.join(sResDir, \
                       "%s_%s.seg%s%d" % (hemi, scROI, segType, segn))
        if bpt2:
            roiResDir += "_pt2"
        
        check_dir(roiResDir)

        fdtp = os.path.join(roiResDir, "fdt_paths.nii.gz")
        check_file(fdtp)

        fdtpn = os.path.join(roiResDir, "fdt_paths_norm.nii.gz")
        check_file(fdtpn)

        print("Processing seed ROI %s, segmentation #%d:" % (scROI, segn))
        for (i1, targROI) in enumerate(h_rois):
            # Get size of ROI
            roi_size_cmd = "fslstats %s -V" % diff_roi_masks[i1]
            # print(roi_size_cmd) # DEBUG
            (stdout, stderr) = Popen(roi_size_cmd.split(' '), \
                                     stdout=PIPE, stderr=PIPE).communicate()
            roiNVoxes = int(stdout.split(" ")[0])
            # print("ROI NVox = %d\n" % roiNVoxes) # DEBUG

            if roiNVoxes == 0:
                print("WARNING: empty ROI mask: %s" % diff_roi_masks[i1])
                continue

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
    resMatFN = os.path.join(sResDir, \
                            "connmats.%s.%s.seg%s.mat" % (hemi, scROI, segType))
    if bpt2:
        resMatFN = resMatFN.replace("connmats.", "connmats.pt2.")    
    if bSpeech:
        resMatFN = resMatFN.replace("connmats.", "connmats.speech.")

    delete_file_if_exists(resMatFN)
    savemat(resMatFN, \
            {"h_rois": h_rois, \
             "connmat_mean": connmat_mean, \
             "connmat_mean_norm": connmat_mean_norm, \
             "connmat_median": connmat_median, \
             "connmat_median_norm": connmat_median_norm})
    check_file(resMatFN)

    print("Connectivity matrices saved to mat file: %s" % resMatFN)

