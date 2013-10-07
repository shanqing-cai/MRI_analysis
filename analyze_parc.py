#!/usr/bin/env python

import os
import sys
import glob
import argparse

import numpy as np

from scai_utils import check_dir, info_log, error_log, check_env_var
from get_qdec_info import get_qdec_info
from parc_stats import get_parc_stats


if __name__ == "__main__":
    ap = argparse.ArgumentParser("Morphological analysis of ROIs based on the specified cortical parcellatio")
    ap.add_argument("fsDir", \
                    help="FreeSurfer subjects dir (e.g., ~/STUT/FSDATA)")
    ap.add_argument("sNameWC", \
                    help="Subject ID wild card (e.g., S??)")
    ap.add_argument("parcName", \
                    help="Parcellation name (e.g., aparc12")
    ap.add_argument("--redo", dest="bRedo", action="store_true", \
                    help="Force redo time-consuming steps")
    ap.add_argument("--skipSubj", \
                    help="List of subjects to skip, separated by commas (e.g., S24,S40")
    #ap.add_argument("groupMethod", \
    #                help="Method for getting the group identity of the subjects (e.g., get_qdec_info,diagnosis). Must be a function name followed by a field name. The field name will be the second input argument to the function. The first input argument will be the subject ID")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    
    #=== Check the matching of environmental SUBJECTS_DIR ===#
    check_env_var("SUBJECTS_DIR", args.fsDir)

    #=== Discover subjects ===#
    check_dir(args.fsDir)
    ds = glob.glob(os.path.join(args.fsDir, args.sNameWC))
    ds.sort()

    skipIDs = []
    if args.skipSubj != None:
        skipIDs = args.skipSubj.split(",")

    sIDs = []
    grps = []
    for (i0, d) in enumerate(ds):
        t_sID = os.path.split(d)[1]

        if skipIDs.count(t_sID) > 0:
            info_log("Skipping subject: %s" % t_sID)
            continue

        sIDs.append(t_sID)
        grps.append(get_qdec_info(t_sID, "diagnosis"))

    ugrps = list((np.unique(np.array(grps))))
    ugrps.sort()

    info_log("Discovered %s subjects" % len(sIDs))
    info_log("The subjects belong to %d groups:" % (len(ugrps)))
    for (i0, grp) in enumerate(ugrps):
        info_log("\t%s" % grp)
                 
    matFile = __file__.replace(".py", ".mat")
    from scipy.io import savemat, loadmat
    
    if not os.path.isfile(matFile) or args.bRedo:
        #=== Extract ROI morphological info ===#
        morphInfo = []
        uniqueROIs = []
        
        for (i0, sID) in enumerate(sIDs):
            morphInfo.append(get_parc_stats(args.fsDir, sID, args.parcName, 
                                            bVerbose=False))
            t_ROIs = list(np.unique(np.array(morphInfo[-1]["rois"])))
            uniqueROIs += t_ROIs

        uniqueROIs = list(np.unique(np.array(uniqueROIs)))

        #=== Construct data structure for analysis ===#
        area_mm2 = {}
        nROIs = len(uniqueROIs)
        for (i0, grp) in enumerate(grps):
            t_idx = np.nonzero(np.array(grps) == grp)[0]
            t_sIDs = list(np.array(sIDs)[t_idx])

            area_mm2[grp] = np.zeros([len(t_sIDs), nROIs])
            
            for (i1, sID) in enumerate(t_sIDs):
                s_idx = sIDs.index(sID)
                t_morphInfo = morphInfo[s_idx]
                t_rois = [x.strip() for x in list(t_morphInfo["rois"])]
            
                for (i2, roi) in enumerate(uniqueROIs):
                    if t_rois.count(roi) == 1:
                        r_idx = t_rois.index(roi)
                        area_mm2[grp][i1, i2] = t_morphInfo["area_mm2"][r_idx]

        #=== Construct and save intermediate results ===#
        intRes = {"sIDs": sIDs, 
                  "grps": grps, 
                  "morphInfo": morphInfo, 
                  "uniqueROIs": uniqueROIs, 
                  "area_mm2": area_mm2}
        savemat(matFile, intRes)
        info_log("Saved intermediate results to file: %s" % matFile)


    info_log("Loaded intermediate results from file: %s" % matFile)
    intRes = loadmat(matFile)

    #== Verify the intermediate results ==#
    bVerified = True;
    bVerified = bVerified and (list(intRes["sIDs"]) == sIDs)
    bVerified = bVerified and (list(intRes["grps"]) == grps)
        
    if bVerified:
        info_log("Intermediate results verified.")
    else:
        error_log("Intermediate results veroficiation failed. Use --redo option to generate valid intermediate results.")

    morphInfo = intRes["morphInfo"]
    uniqueROIs = list(intRes["uniqueROIs"])

    uniqueROIs = [x.strip() for x in uniqueROIs]

    area_mm2 = intRes["area_mm2"][0]
    
    #=== Test analysis ===#
    import matplotlib.pyplot as plt

    roi_a = "rh_vPMC"
    roi_b = "rh_PT"
    
    ridx_a = uniqueROIs.index(roi_a)
    ridx_b = uniqueROIs.index(roi_b)

    x_PFS_a = area_mm2["PFS"][0][:, ridx_a]
    x_PFS_b = area_mm2["PFS"][0][:, ridx_b]
    x_PWS_a = area_mm2["PWS"][0][:, ridx_a]
    x_PWS_b = area_mm2["PWS"][0][:, ridx_b]

    plt.figure()
    plt.plot(x_PFS_a, x_PFS_b, "bo")
    plt.plot(x_PWS_a, x_PWS_b, "ro")
    plt.xlabel("Surface area of %s" % roi_a)
    plt.ylabel("Surface area of %s" % roi_b)
    plt.show()
