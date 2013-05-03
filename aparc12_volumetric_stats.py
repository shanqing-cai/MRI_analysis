#!/usr/bin/python

import os
import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

from get_qdec_info import get_qdec_info
from aparc12 import *
from scai_utils import *

BASE_DIR = "/users/cais/STUT/analysis/aparc12_tracts"
DATA_DIR = "/users/cais/STUT/DATA"
CTAB = "/users/cais/STUT/slaparc_550.ctab"

SEGSTATS_SUM_WC = "aparc12_wm%dmm.segstats.txt"

hemis = ["lh", "rh"]
grps = ["PFS", "PWS"]
grpColors = {"PFS": [0, 0, 0], "PWS": [1, 0, 0]}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Analyze volumes of GM and WM volumes in aparc12")
    ap.add_argument("wmDepth", type=int, help="WM depth in mm (e.g., 1, 2)")
    ap.add_argument("matter", type=str, choices=["GM", "WM"], help="Gray or white matter")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    # === Args input arguments === 
    args = ap.parse_args()
    wmDepth = args.wmDepth
    matter = args.matter
        
    # === Determine the subject list and their group memberships ===
    check_dir(BASE_DIR)
    ds = glob.glob(os.path.join(BASE_DIR, "S??"))
    ds.sort()

    sIDs = []
    isPWS = []
    for (i0, t_path) in enumerate(ds):
        (t_path_0, t_sID) = os.path.split(t_path)
        sIDs.append(t_sID)

        if get_qdec_info(t_sID, "diagnosis") == "PWS":
            isPWS.append(1)
        else:
            isPWS.append(0)
        
    isPWS = np.array(isPWS)

    assert(len(sIDs) > 0)
    assert(len(sIDs) == len(isPWS))

    # === Get the list of cortical ROIs ===
    rois0 = get_aparc12_cort_rois(bSpeech=True)

    check_file(CTAB)
    (ctab_roi_nums, ctab_roi_names) = read_ctab(CTAB)

    # Duplex into both hemispheres
    roi_names = []
    roi_nums = []

    for (i0, hemi) in enumerate(hemis):
        for (i1, roi) in enumerate(rois0):
            t_roi_name = "%s_%s" % (hemi, roi)
            
            if matter == "WM":
                t_roi_name += "_wm"

            assert(ctab_roi_names.count(t_roi_name) == 1)
            
            idx = ctab_roi_names.index(t_roi_name)
            roi_names.append(t_roi_name)
            roi_nums.append(ctab_roi_nums[idx])
        
    
    # === Loop through all subjects === #
    assert(len(roi_names) == len(roi_nums))

    nROIs = len(roi_names)
    ns = len(sIDs)

    volGM = np.zeros([ns, nROIs])

    for (i0, sID) in enumerate(sIDs):
        sDataDir = os.path.join(DATA_DIR, sID)
        check_dir(sDataDir)

        sumfn = os.path.join(sDataDir, SEGSTATS_SUM_WC % wmDepth)
        check_file(sumfn)

        sumt = remove_empty_strings(read_text_file(sumfn))
        
        t_labs = []
        t_vol_mm3 = []
        for (i1, tline) in enumerate(sumt):
            titems = remove_empty_strings(tline.replace('\t', ' ').split(' '))
            
            if titems[0] == "#":
                continue
            assert(len(titems) == 5)
                
            t_labs.append(int(titems[1]))
            t_vol_mm3.append(float(titems[3]))

        for (i1, t_roi_num) in enumerate(roi_nums):
            if t_labs.count(t_roi_num) == 1:
                volGM[i0, i1] = t_vol_mm3[t_labs.index(t_roi_num)]
            else:
                print("WARNING: label %s missing in subject %s's file: %s" % \
                      (t_roi_num, sID, sumfn))


    # === Statistical comparison === #
    mean_volGM = {}
    ste_volGM = {}
    nsg = {}

    nsg["PFS"] = len(np.nonzero(isPWS == 0))
    nsg["PWS"] = len(np.nonzero(isPWS == 1))

    mean_volGM["PFS"] = np.mean(volGM[isPWS == 0], axis=0)
    ste_volGM["PFS"] = np.std(volGM[isPWS == 0], axis=0) / np.sqrt(nsg["PFS"])
    
    mean_volGM["PWS"] = np.mean(volGM[isPWS == 1], axis=0)
    ste_volGM["PWS"] = np.std(volGM[isPWS == 1], axis=0) / np.sqrt(nsg["PWS"])

    # === Visualiation === #
    for (i0, grp) in enumerate(grps):
        plt.errorbar(range(nROIs), mean_volGM[grp], yerr=ste_volGM[grp], \
                     color=grpColors[grp])
    plt.xticks(range(nROIs), roi_names, rotation=90.0)
    plt.show()
