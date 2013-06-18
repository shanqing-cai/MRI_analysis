#!/usr/bin/python

import os
import sys
import glob
import argparse

from aparc12 import get_aparc12_cort_rois
from scai_utils import *

tractSegDir = "/users/cais/STUT/analysis/tractseg_aparc12/"
L2_DIR = "/users/cais/STUT/analysis/tract_subcort_conn"

allROIs = get_aparc12_cort_rois(lobe="all", bSpeech=False)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Batch run of tractography- and aparc12-based subcorticla-cortical connectivity analysis, by calling tract_seg.py and tract_subcort_conn.py")
    ap.add_argument("scSeed", \
                    help="Subcortical seed (e.g., lh.Thalamus-Proper)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    check_dir(tractSegDir)

    ds = glob.glob(os.path.join(tractSegDir, "S??"))
    ds.sort()

    sIDs = []

    for (i0, t_fn) in enumerate(ds):
        (t_path, t_sID) = os.path.split(t_fn)
        sIDs.append(t_sID)

    # === Check the existence of the subcortical seed mask file === %
    scMask = os.path.join(L2_DIR, "%s_ccStop" % args.scSeed, \
                          "mask", "merged_mask_min_bin.nii.gz")
    check_file(scMask)
        
    for roi in allROIs:
        for sID in sIDs:
            saydo("tract_seg.py %s %s roiconn --roi %s --ccStop" % \
                  (sID, args.scSeed, roi))

        saydo("tract_subcort_conn.py %s --ccstop --roi %s --mask %s" % \
              (args.scSeed, roi, scMask))


    
    
