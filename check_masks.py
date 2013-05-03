#!/usr/bin/python

import os
import sys
import argparse
import tempfile

from scai_utils import *

CTAB = "/users/cais/STUT/slaparc_550.ctab"
DATA_DIR = "/users/cais/STUT/DATA"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Manual checking of structural masks of SLaparc ROIs. Depends on the files such as aparc12_wm1mm.nii.gz")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("wmDepth", type=int, help="WM depth (e.g., 1, 2)")
    ap.add_argument("regionNames", help="Names of the regions, separated by commas (e.g., lh_pSTg_wm,lh_pdSTs_wm)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    # === Parse input arguments ===
    args = ap.parse_args()
    sID = args.sID
    wmDepth = args.wmDepth
    regionNames = args.regionNames

    # === Locate files === 
    sDataDir = os.path.join(DATA_DIR, sID)
    check_dir(sDataDir)
    
    vol = os.path.join(sDataDir, "aparc12_wm%dmm.nii.gz" % wmDepth)
    check_file(vol)

    ROIs = regionNames.split(',')
    ROINums = [-1] * len(ROIs)

    # === Read color table (ctab) ===
    check_file(CTAB)
    (roi_nums, roi_names) = read_ctab(CTAB)
    
    mask_fns = [""] * len(ROIs)
    for (i0, t_roi) in enumerate(ROIs):
        assert(roi_names.count(t_roi) == 1)
        idx = roi_names.index(t_roi)
        ROINums[i0] = roi_nums[idx]

        mask_fns[i0] = tempfile.mktemp() + "_%s.nii.gz" % t_roi
        bin_cmd = "mri_binarize --i %s --match %d --o %s" % \
                  (vol, roi_nums[idx], mask_fns[i0])
        saydo(bin_cmd)
        check_file(mask_fns[i0])        

    # === Let user view the images === #
    view_cmd = "freeview %s " % vol
    for (i0, t_maskfn) in enumerate(mask_fns):
        view_cmd += "%s:colormap=jet " % t_maskfn
    
    saydo(view_cmd)
    
    for (i0, t_maskfn) in enumerate(mask_fns):
        saydo("rm -f %s" % t_maskfn)
    
