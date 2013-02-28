#!/usr/bin/python

import os
import sys
import glob
import argparse
import tempfile

from aparc12 import get_aparc12_cort_rois
from scai_utils import *
from aparc12_probtrackx_2 import get_roi_num

#TRACTS_RES_DIR = "/users/cais/STUT/analysis/aparc12_tracts_2"
FSDATA_DIR = "/users/cais/STUT/FSDATA"
TRACULA_BASE = "/users/cais/STUT/analysis/dti2/tracula"

CTAB_FN_CORTICAL = '/users/cais/STUT/slFRS17.ctab'


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate the diffusion-space gray-matter mask of an ROI in aparc12 of a subject")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("hemi", help="Hemisphere {lh, rh}")
    ap.add_argument("roi", help="Name of the aparc12 cortical ROI")
    ap.add_argument("outFN", help="File name of the output diffusion-space GM mask")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
        
    args = ap.parse_args()
    sID = args.sID
    hemi = args.hemi
    roi = args.roi
    outFN = args.outFN

    # Locate the aparc12 volume
    a12v = os.path.join(FSDATA_DIR, sID, "mri", "aparc12.nii.gz")
    check_file(a12v)

    # Get the number of the ROI
    roiNum = get_roi_num(CTAB_FN_CORTICAL, hemi, roi)

    # Locate the anat-to-diff xfm file
    sTracDir = os.path.join(TRACULA_BASE, sID)
    check_dir(sTracDir)

    anat2diff_xfm_mat = os.path.join(sTracDir, "dmri", "xfms", \
                                     "anatorig2diff.bbr.mat") 
    check_file(anat2diff_xfm_mat)

    # Locate the FA img
    FA_img = os.path.join(sTracDir, "dmri", "dtifit_FA.nii.gz")
    check_file(FA_img)
    
    # Binarize
    binOut = tempfile.mktemp() + ".nii.gz"
    binCmd = "mri_binarize --i %s --match %d --o %s" % \
             (a12v, roiNum, binOut)
    saydo(binCmd)
    check_file(binOut)
    
    # Transform
    xfm_cmd = "flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour"\
               %(binOut, FA_img, anat2diff_xfm_mat, outFN)
    saydo(xfm_cmd)
    check_file(outFN)
       
    # Clean up temporary files
    os.system("rm -f %s" % binOut)

    # Show command line for checking registration
    checkCmd = "tkregister2 --mov %s --targ %s --identity --reg identity.dat" \
               % (outFN, FA_img)
    print("Command line for checking the output mask file: %s\n" % outFN)
    print("\t%s" % checkCmd)
