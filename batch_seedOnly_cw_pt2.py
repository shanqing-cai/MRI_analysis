#!/usr/bin/python

import os
import sys
import argparse

from scai_utils import *

"""
speechROIs = ['H', 'PO', 'PP', 'PT', 'SPL', 'TP',  'SMA', 'aCG', 'aCO', \
              'aFO', 'aINS', 'aSTg', 'adSTs', \
              'aSMg', 'dIFo', 'dMC', 'dSC', 'mdPMC', 'midMC', 'midPMC', 'pCO', \
              'pFO', 'pIFs', 'pINS', 'aSTg', 'pSTg', 'pdPMC', 'pdSTs', \
              'preSMA', 'vIFo', 'vMC', 'vPMC', 'vSC', 'midCG']
"""

from aparc12_oldVer import get_aparc12_cort_rois
speechROIs = get_aparc12_cort_rois("all", bSpeech=True)

#speechROIs = ['midCG']

BIN = "/users/cais/STUT/scripts/aparc12_probtrackx_2.py"
RES_DIR = "/users/cais/STUT/analysis/aparc12_tracts_pt2"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run aparc12_probtrackx_2.py under seedOnly, caww (corpus-callosum avoidance, ipsilateral WM waypoint) mode on all speech-network ROIs of a hemisphere of a subject")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("hemi", help="Hemisphere (e.g., lh or lh,rh)")
    ap.add_argument("--noSurfProj", dest="bNoSurfProj", action="store_true", \
                    help="Do not perform projection to surface")
    ap.add_argument("--redo", dest="bRedo", action="store_true", \
                    help="Force redo all ROIs. If False, then the script will skip the ROIs that have already been completed.")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    sID = args.sID
    hemi = args.hemi
    bRedo = args.bRedo

    if args.hemi.count(",") == 0:
        hemis = [hemi]
    else:
        hemis = args.hemi.split(",")

    sDir = os.path.join(RES_DIR, sID)
    
    for (j0, t_hemi) in enumerate(hemis):
        for (i0, roi) in enumerate(speechROIs):
            # Check if the current ROI is already done
            roiDir = os.path.join(sDir, "%s_%s_cw_pt2" % (t_hemi, roi))
            bDone = os.path.isfile(os.path.join(roiDir, "fdt_paths.nii.gz")) \
                    and os.path.isfile(os.path.join(roiDir, "seed_size.txt")) \
                    and os.path.isfile(os.path.join(roiDir, "waytotal")) \
                    and os.path.isfile(os.path.join(roiDir, "viewcmds.txt"))
            if not args.bNoSurfProj:
                bDone = bDone and (os.path.isfile(os.path.join(roiDir, "fdt_paths_norm.lh.mgz")) \
                                   and os.path.isfile(os.path.join(roiDir, "fdt_paths_norm.rh.mgz")))

            if bDone and (not bRedo):
                print("INFO: It appears that ROI %s has already been processed. Skipping it because bRedo = False." % (roi))

            if (not bDone) or (bRedo):
                cmd = "%s %s %s_%s --pt2 --ccwaypoint" % \
                      (BIN, sID, t_hemi, roi)

                if args.bNoSurfProj:
                    cmd += " --noSurfProj"

                saydo(cmd)
