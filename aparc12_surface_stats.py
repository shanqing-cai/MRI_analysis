#!/usr/bin/python

import os
import sys
import glob
import argparse
import tempfile
import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.stats as stats
from copy import deepcopy
from subprocess import Popen, PIPE

from get_qdec_info import get_qdec_info
from fs_load_stats import fs_load_stats
from aparc12 import *
from scai_utils import *
from scai_stats import cohens_d

BASE_DIR = "/users/cais/STUT/analysis/aparc12_tracts"
DATA_DIR = "/users/cais/STUT/DATA"
FSDATA_DIR = "/users/cais/STUT/FSDATA"
CTAB = "/users/cais/STUT/slaparc_550.ctab"

SEGSTATS_SUM_WC = "aparc12_wm%dmm.segstats.txt"

P_THRESH_UNC = 0.05

hemis = ["lh", "rh"]
grps = ["PFS", "PWS"]
grpColors = {"PFS": [0, 0, 0], "PWS": [1, 0, 0]}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Analyze aparc12 surface annotation: Surface area and average thickness")
    ap.add_argument("-r", dest="bReload", action="store_true", \
                    help="Reload data (time-consuming)")
    # ap.add_argument("hemi", help="Hemisphere {lh, rh}")
    
    # if len(sys.argv) == 1:
    #     ap.print_help()
    #     sys.exit(0)

    # === Args input arguments === #
    args = ap.parse_args()
    bReload = args.bReload
    # hemi = args.hemi
    # assert(hemis.count(hemi) == 1)

    # === Determine the subject list and their group memberships === #
    check_dir(BASE_DIR)
    ds = glob.glob(os.path.join(BASE_DIR, "S??"))
    ds.sort()

    sIDs = []
    isPWS = []
    SSI4 = []
    for (i0, t_path) in enumerate(ds):
        (t_path_0, t_sID) = os.path.split(t_path)
        sIDs.append(t_sID)
        SSI4.append(get_qdec_info(t_sID, "SSI"))

        if get_qdec_info(t_sID, "diagnosis") == "PWS":
            isPWS.append(1)
        else:
            isPWS.append(0)
        
    isPWS = np.array(isPWS)
    SSI4 = np.array(SSI4)

    assert(len(sIDs) > 0)
    assert(len(sIDs) == len(isPWS))

    # === Get the list of cortical ROIs (Speech network only) ===
    rois0 = get_aparc12_cort_rois(bSpeech=True)

    check_file(CTAB)
    (ctab_roi_nums, ctab_roi_names) = read_ctab(CTAB)

    # Duplex into both hemispheres
    roi_names = []
    roi_nums = []

    for (i0, hemi) in enumerate(hemis):
        for (i1, roi) in enumerate(rois0):
            t_roi_name = "%s_%s" % (hemi, roi)
            
            assert(ctab_roi_names.count(t_roi_name) == 1)
            
            idx = ctab_roi_names.index(t_roi_name)
            roi_names.append(t_roi_name)
            roi_nums.append(ctab_roi_nums[idx])

    assert(len(roi_names) == len(roi_nums))
    
    # === Load data: Loop through all subjects === #
    cachePklFN = "aparc12_surface_stats_dset.pkl"

    nROIs = len(roi_names)
    ns = len(sIDs)


    if bReload:
        print("INFO: bReload = True: Reloading data (time-consuming)\n")

        labArea = np.zeros([ns, nROIs])
        labAvgThick = np.zeros([ns, nROIs])
        # Label area normalized by hemisphere surface area
        labAreaNorm = np.zeros([ns, nROIs])

        for (i0, sID) in enumerate(sIDs):
            t_rois = []
            t_roi_nums = []
            t_area = []
            t_area_norm = []
            t_thick = []

            for (i1, hemi) in enumerate(hemis):
                # == Load hemisphere total surface area == #
                hemiStatsFN = os.path.join(FSDATA_DIR, sID, \
                                           "stats", "%s.aparc.stats" % hemi)
                check_file(hemiStatsFN)
                t_hemiSurfArea = fs_load_stats(hemiStatsFN, "SurfArea")

                tmpfn = tempfile.mktemp()

                hthick = os.path.join(FSDATA_DIR, sID, "surf", \
                                      "%s.thickness" % hemi)
                check_file(hthick)

                print("Loading data from subject %s, hemisphere %s" \
                      % (sID, hemi))
                    
                (sout, serr) = Popen(["mri_segstats", "--annot", \
                                      sID, hemi, "aparc12", \
                                      "--in", hthick, \
                                      "--sum", tmpfn], \
                                     stdout=PIPE, stderr=PIPE).communicate()

                sout = read_text_file(tmpfn)
                os.system("rm -rf %s" % tmpfn)

                k0 = 0
                while (sout[k0].startswith("# ")):
                    k0 = k0 + 1
            
                sout = sout[k0 :]

                for tline in sout:
                    if len(tline) == 0:
                        break;
                
                    t_items = remove_empty_strings(\
                        tline.replace('\t', ' ').split(' '))

                    if len(t_items) == 10:
                        t_rois.append(hemi + "_" + t_items[4])
                        
                        if hemi == "lh":
                            t_roi_nums.append(1000 + int(t_items[1]))
                        else:
                            t_roi_nums.append(2000 + int(t_items[1]))

                        t_area.append(float(t_items[3]))
                        t_area_norm.append(float(t_items[3]) / t_hemiSurfArea)
                        t_thick.append(float(t_items[5]))

            # == Matching and filling values == #
            for (i2, t_rn) in enumerate(roi_nums):
                if t_roi_nums.count(t_rn) > 0:
                    idx = t_roi_nums.index(t_rn)
                    labArea[i0][i2] = t_area[idx]
                    labAreaNorm[i0][i2] = t_area_norm[idx]
                    labAvgThick[i0][i2] = t_thick[idx]

        # === Save to pkl file === #        
        dset = {"labArea": labArea, \
                "labAreaNorm": labAreaNorm, \
                "labAvgThick": labAvgThick}

        os.system("rm -rf %s" % cachePklFN)
        cachePklF = open(cachePklFN, "wb")
        pickle.dump(dset, cachePklF)
        cachePklF.close()
        check_file(cachePklFN)

        print("INFO: Saved loaded data to file: %s\n" % os.path.abspath(cachePklFN))
    else:
        print("INFO: Loading saved data from file: %s\n" % os.path.abspath(cachePklFN))

        cachePklF = open(cachePklFN, "rb")
        dset = pickle.load(cachePklF)
        cachePklF.close()

        labArea = dset["labArea"]
        labAreaNorm = dset["labAreaNorm"]
        labAvgThick = dset["labAvgThick"]

    # === Check data validity === #
    assert(len(labArea) == ns)
    assert(len(labAreaNorm) == ns)
    assert(len(labAvgThick) == ns)

    # === Statistical comparison === #
    mean_area = {}
    std_area = {}
    nsg = {}

    for (i0, grp) in enumerate(grps):
        nsg[grp] = len(np.nonzero(isPWS == i0))

        mean_area[grp] = np.mean(labArea[isPWS == i0], axis=0)
        # std_area[grp] = np.std(labArea[isPWS == i0], axis=0) / np.sqrt(nsg[grp])
        std_area[grp] = np.std(labArea[isPWS == i0], axis=0)

    cmprItems = ["labArea", "labAreaNorm", "labAvgThick"]
    for (h0, cmprItem) in enumerate(cmprItems):
        print("--- List of significant differences in %s (p < %f uncorrected) ---"  \
              % (cmprItem, P_THRESH_UNC))
        p_tt_val = np.zeros([nROIs])
        t_tt_val = np.zeros([nROIs])

        for (i0, t_roi) in enumerate(roi_names):
            if h0 == 0:
                dat_PWS = labArea[isPWS == 1, i0]
                dat_PFS = labArea[isPWS == 0, i0]
            elif h0 == 1:
                dat_PWS = labAreaNorm[isPWS == 1, i0]
                dat_PFS = labAreaNorm[isPWS == 0, i0]
            elif h0 == 2:
                dat_PWS = labAvgThick[isPWS == 1, i0]
                dat_PFS = labAvgThick[isPWS == 0, i0]
                
            (t_tt, p_tt) = stats.ttest_ind(dat_PWS, dat_PFS)
        
            p_tt_val[i0] = p_tt
            t_tt_val[i0] = t_tt

            if p_tt_val[i0] < P_THRESH_UNC:
                if t_tt_val[i0] < 0:
                    dirString = "PWS < PFS"
                else:
                    dirString = "PWS > PFS"            
                print("%s: p = %f; t = %f (%s)" \
                      % (t_roi, p_tt_val[i0], t_tt_val[i0], dirString))

                print("\tMean +/- SD: PWS: %.5f +/- %.5f; PFS: %.5f +/- %.5f" \
                      % (np.mean(dat_PWS), np.std(dat_PWS), \
                         np.mean(dat_PFS), np.std(dat_PFS)))
                print("\tCohens_d = %.3f" % cohens_d(dat_PWS, dat_PFS))
                
        print("\n")
    

    # === Spearman correlation === #
    for (h0, cmprItem) in enumerate(cmprItems):
        print("--- Spearman correlations with SSI4 in %s (p < %f uncorrected) ---"  \
              % (cmprItem, P_THRESH_UNC))
        p_spc_val = np.zeros([nROIs])
        rho_spc_val = np.zeros([nROIs])

        for (i0, t_roi) in enumerate(roi_names):
            if h0 == 0:
                (r_spc, p_spc) = stats.spearmanr(SSI4[isPWS == 1], \
                                               labArea[isPWS == 1, i0])
            elif h0 == 1:
                (r_spc, p_spc) = stats.spearmanr(SSI4[isPWS == 1], \
                                                labAreaNorm[isPWS == 1, i0])
            elif h0 == 2:
                (r_spc, p_spc) = stats.spearmanr(SSI4[isPWS == 1], \
                                                labAvgThick[isPWS == 1, i0])
        
            p_spc_val[i0] = p_spc
            rho_spc_val[i0] = r_spc

            if p_spc_val[i0] < P_THRESH_UNC:
                if rho_spc_val[i0] < 0:
                    dirString = "-"
                else:
                    dirString = "+"

                print("%s: p = %f; rho = %f (%s)" \
                      % (t_roi, p_spc_val[i0], rho_spc_val[i0], dirString))
                
        print("\n")

    # === Compare combined dIFo and vIFo === #
    lh_IFo_area = {}
    lh_IFo_areaNorm = {}
    for (i0, grp) in enumerate(grps):
        lh_IFo_area[grp] = labArea[isPWS == i0, roi_names.index("lh_vIFo")] + \
                           labArea[isPWS == i0, roi_names.index("lh_dIFo")]
        lh_IFo_areaNorm[grp] = labAreaNorm[isPWS == i0, roi_names.index("lh_vIFo")] + \
                               labAreaNorm[isPWS == i0, roi_names.index("lh_dIFo")]

    (t_tt, p_tt) = stats.ttest_ind(lh_IFo_area["PWS"], \
                                   lh_IFo_area["PFS"])
    print("-- Comparing lh_IFo area: --")
    print("\tp = %f; t = %f" % (p_tt, t_tt))
    print("\tPWS: %.1f +/- %.1f; PFS: %.1f +/- %.1f" \
          % (np.mean(lh_IFo_area["PWS"]), np.std(lh_IFo_area["PWS"]), \
             np.mean(lh_IFo_area["PFS"]), np.std(lh_IFo_area["PFS"])))
    print("\n")

    (t_tt, p_tt) = stats.ttest_ind(lh_IFo_areaNorm["PWS"], \
                                   lh_IFo_areaNorm["PFS"])
    print("-- Comparing lh_IFo areaNorm: --")
    print("\tp = %f; t = %f" % (p_tt, t_tt))
    print("\tPWS: %.1e +/- %.1e; PFS: %.1e +/- %.1e" \
          % (np.mean(lh_IFo_areaNorm["PWS"]), np.std(lh_IFo_areaNorm["PWS"]), \
             np.mean(lh_IFo_areaNorm["PFS"]), np.std(lh_IFo_areaNorm["PFS"])))
    

    # === Correlating combined IFo with SSI4 === #
    (r_spc, p_spc) = stats.spearmanr(SSI4[isPWS == 1], lh_IFo_area["PWS"])
    print("-- Correlating SSI4 with lh_IFo area: --")
    print("\tp = %f; rho = %f" % (p_spc, r_spc))
    print("\n")

    (r_spc, p_spc) = stats.spearmanr(SSI4[isPWS == 1], lh_IFo_areaNorm["PWS"])
    print("-- Correlating SSI4 with lh_IFo areaNorm: --")
    print("\tp = %f; rho = %f" % (p_spc, r_spc))
    print("\n")
    
        
    # === Visualiation === #
    """
    for (i0, grp) in enumerate(grps):
       plt.errorbar(range(nROIs), mean_area[grp], yerr=std_area[grp], \
                    color=grpColors[grp])
    plt.xticks(range(nROIs), roi_names, rotation=90.0)
    plt.show()
    """


