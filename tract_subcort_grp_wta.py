#!/usr/bin/python

import os
import sys
import glob
import argparse

from scai_utils import *

TRACT_SUBCORT_DIR = "/users/cais/STUT/analysis/tract_subcort_conn"
TRACT_SUBCORT_WTA_DIR = "/users/cais/STUT/analysis/tract_subcort_wta"

MODES = ["coarse", "fine"]

REGIONS = ['Prefrontal', 'Premotor', 'Precentral', 'Postcentral', \
           'PPC', 'Occipital', 'Temporal']
GRPS = ["PWS", "PFS"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Group-level mean-based WTA parcellation of a subcortical nucleus (requires mask)")
    ap.add_argument("scSeed", help="Name of the subcortical nucleus (e.g., lh.Putamen_ccStop")
    ap.add_argument("mode", help="Mode {coarse | fine}")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    # Parse input arguments
    args = ap.parse_args()
    assert(MODES.count(args.mode) == 1)
    
    check_dir(TRACT_SUBCORT_DIR)
    scDir = os.path.join(TRACT_SUBCORT_DIR, args.scSeed)
    check_dir(scDir)

    check_dir(TRACT_SUBCORT_WTA_DIR, bCreate=True)
    outDir = os.path.join(TRACT_SUBCORT_WTA_DIR, args.mode)
    check_dir(outDir, bCreate=True)
    
    # Look for the mask
    mask = os.path.join(TRACT_SUBCORT_DIR, args.scSeed, \
                        "mask", "merged_mask_min_bin.nii.gz")
    check_file(mask)
    print("INFO: Mask = %s" % mask)
    
    # Collect and merge the data 
    ROIs = []
    if args.mode == "coarse":
        for (i0, t_region) in enumerate(REGIONS):
            ROIs.append("%d_%s" % (i0, t_region))

    meanMerged = {}
    for (i0, grp) in enumerate(GRPS):
        meanMerged[grp] = os.path.join(outDir, "meanMerged_%s.nii.gz" % grp)
        mergeCmd = "fslmerge -t %s " % meanMerged[grp]
        for (i1, roi) in enumerate(ROIs):    
            grpMergedWC = os.path.join(TRACT_SUBCORT_DIR, args.scSeed, roi, \
                                       "merged_*_mean_%s.nii.gz" % grp)
            grpMerged = glob.glob(grpMergedWC)
            assert(len(grpMerged) == 1)
            grpMerged = grpMerged[0]
            
            mergeCmd += "%s " % grpMerged

        saydo(mergeCmd)
        check_file(meanMerged[grp])

        maxnFN = os.path.join(outDir, "maxn_%s.nii.gz" % grp)
        maxnCmd = "fslmaths %s -Tmaxn %s" % (meanMerged[grp], maxnFN)
        saydo(maxnCmd)
        check_file(maxnFN)

        (sout, serr) = cmd_stdout("fslstats %s -P 100" % mask)
        assert(len(serr) == 0)
        assert(float(sout.split(" ")[0]) == 1.0)r

        addCmd = "fslmaths %s -add %s %s" % (maxnFN, mask, maxnFN)
        saydo(addCmd)
        
    

        
    
    
