#!/usr/bin/python

import os
import sys
import glob
import argparse

from get_qdec_info import get_qdec_info
from scai_utils import *

DATA_DIR = "/users/cais/STUT/DATA"
TRACULA_DIR = "/users/cais/STUT/analysis/dti2/tracula"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate probabilistic map in TBSS space for aparc12 (SLaparc) parcellations")
    ap.add_argument("tbssDir", \
                    help="TBSS directory (e.g., /users/cais/STUT/analysis/dt_tbss_dtiprep2)")
    ap.add_argument("--PFSOnly", dest="bPFSOnly", action="store_true", \
                    help="Use only PFS (control) subjects to generate the probabilistic map")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
    
    args = ap.parse_args()
    tbssDir = args.tbssDir
    bPFSOnly = args.bPFSOnly

    # === Get subject list === #
    check_dir(tbssDir)
    origDir = os.path.join(tbssDir, "origdata")
    check_dir(origDir)

    aparc12Dir = os.path.join(tbssDir, "aparc12")
    check_dir(aparc12Dir, bCreate=True)

    ds = glob.glob(os.path.join(origDir, "S??.nii.gz"))
    ds.sort()
    sIDs = []
    isPWS = []
    aparc12_fns = []
    aparc12_diff_fns = []

    merged = os.path.join(aparc12Dir, "merged.nii.gz")
    os.system("rm -f %s" % merged)
    merge_cmd = "fslmerge -t %s " % merged
    
    for (i0, d) in enumerate(ds):
        [tpath, tfn] = os.path.split(d)
        sID = tfn.replace(".nii.gz", "")
        sIDs.append(sID)
        
        isPWS.append(get_qdec_info(sID, "diagnosis") == "PWS")
        
        t_aparc12 = os.path.join(DATA_DIR, sID, "aparc12.nii.gz")
        check_file(t_aparc12)
        aparc12_fns.append(t_aparc12)

        # == Locate the d2a FSL xfm mat == #
        d2a = os.path.join(TRACULA_DIR, sID, "dmri", "xfms", "d2a.mat")
        check_file(d2a)
        
        # == Use convert_xfm to create a2d FSL xfm mat == #
        a2d = os.path.join(TRACULA_DIR, sID, "dmri", "xfms", "a2d.mat")
        os.system("rm -f %s" % a2d)
        inv_cmd = "convert_xfm -omat %s -inverse %s" % (a2d, d2a)
        saydo(inv_cmd)
        check_file(a2d)
        
        aparc12_diff_fns.append(os.path.join(aparc12Dir, \
                               "%s_aparc12_diff.nii.gz" % sID))

        t_img = os.path.join(aparc12Dir, "%s_tmp.nii.gz" % sID)
        os.system("rm -f %s" % aparc12_diff_fns[-1])
        xfm_cmd = "flirt -in %s -ref %s -applyxfm -init %s -o %s" % \
                 (t_aparc12, d, a2d, t_img) + \
                 " -interp nearestneighbour"
        saydo(xfm_cmd)
        check_file(t_img)

        # == Transform to common diffusion space == #
        premat = os.path.join(tbssDir, "FA", "%s_FA_to_target.mat" % sID)
        check_file(premat)

        warp = os.path.join(tbssDir, "FA", "%s_FA_to_target_warp.nii.gz" % sID)
        check_file(warp)

        ref = os.path.join(tbssDir, "stats", "mean_FA.nii.gz")
        check_file(ref)

        os.system("rm -f %s" % aparc12_diff_fns[-1])
        warpCmd = "applywarp --ref=%s --in=%s " % (ref, t_img) + \
                  "--warp=%s " % (warp) + \
                  "--out=%s --interp=nn" % (aparc12_diff_fns[-1])
        saydo(warpCmd)
        check_file(aparc12_diff_fns[-1])

        merge_cmd += "%s " % aparc12_diff_fns[-1]

    # === Concatenate === #
    saydo(merge_cmd)
    check_file(merged)

    print("INFO: Merged aparc12 label file generated at: %s" % merged)
    

        
    
