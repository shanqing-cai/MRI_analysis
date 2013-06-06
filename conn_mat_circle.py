#!/usr/bin/python

import os
import sys
import argparse
import numpy as np
import pylab as pl
import scipy.io 
from copy import deepcopy
from scai_mne.viz import circular_layout, plot_connectivity_circle

from scai_utils import *
from aparc12 import get_aparc12_cort_rois

lobes = ["Prefrontal", "Premotor", "Insular", "Precentral", \
         "Postcentral", "PPC", "Temporal", "Cingulate"]
# lobeClrs = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 0, 0), \
#            (1, 1, 0), (0, 0.5, 0), (1, 0.5, 0), (0.5, 0, 0.5)]
lobeClrs = [(0.5, 0.5, 0.5)] * len(lobes)


COORD_FILE = "/users/cais/STUT/FSDATA/fsaverage2/mri/aparc12_roi_coords.txt"

hemis=["lh", "rh"]

FIG_DIR = "/users/cais/STUT/figures"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Draw connectivity circle plot")
    ap.add_argument("inMatFN", help="Input mat file with the a_cmat")
    ap.add_argument("hemi", type=str, choices=hemis, help="Hemisphere")
    ap.add_argument("grp", type=str, help="Group (e.g., PWS, PFS: must exist as a_cmat[grp] in inMatFN")
    ap.add_argument("--vmax", type=float, default=np.nan,
                    help="Maximum value (e.g., 331.8")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    # === Parse input arguments === #
    args = ap.parse_args()
    inMatFN = args.inMatFN
    hemi = args.hemi
    grp = args.grp
    vmax = args.vmax

    # === ROIs by lobe ===
    rois_bl = {}    
    for (i0, t_lobe) in enumerate(lobes):
        rois_bl[t_lobe] = get_aparc12_cort_rois(lobe=t_lobe, bSpeech=True)
        rois_bl[t_lobe] = np.array(rois_bl[t_lobe])

    # === Read the ROI centers of gravity from text file === #
    # check_file(COORD_FILE)

    cf = open(COORD_FILE, "rt")
    ct = cf.read().split('\n')
    ct = remove_empty_strings(ct)
    cf.close()

    roi_names = []
    roi_nums = []
    roi_coords = []
    for (i0, tline) in enumerate(ct):
        t_items = tline.split(' ')
        if len(t_items) != 5:
            raise Exception, "Unrecognized formant in a line of %s: %s" \
                             % (COORD_FILE, tline)
        roi_names.append(t_items[0])
        roi_nums.append(t_items[1])
        t_coord = [float(t_items[2]), float(t_items[3]), float(t_items[4])]
        roi_coords.append(t_coord)


    cogy_bl = {}
    for (i0, t_lobe) in enumerate(lobes):
        cogy_bl[t_lobe] = np.zeros(len(rois_bl[t_lobe]))

        for (i1, t_roi) in enumerate(rois_bl[t_lobe]):
            assert(roi_names.count("lh_" + t_roi) == 1)

            t_coord = roi_coords[roi_names.index("lh_" + t_roi)]

            cogy_bl[t_lobe][i1] = t_coord[1]

            # print("%s - %f" % (t_roi, t_coord[1])) # DEBUG
            
        sortidx = sorted(range(len(cogy_bl[t_lobe])), \
                         key=lambda k: cogy_bl[t_lobe][k], \
                         reverse=True)
        rois_bl[t_lobe] = rois_bl[t_lobe][sortidx]

    # === Combine into a single list of ROIs === #
    rois = []
    
    for (i0, t_lobe) in enumerate(lobes):
        rois += rois_bl[t_lobe]

    for (i0, t_roi) in enumerate(rois):
        rois[i0] = hemi[0].upper() + " " + t_roi

    rois = np.array(rois)
    nrois = len(rois)

    roi_clrs = [()] * nrois
    ccnt = 0
    for (i0, t_lobe) in enumerate(lobes):
        for i1 in range(len(rois_bl[t_lobe])):
            roi_clrs[ccnt] = lobeClrs[i0]
            ccnt += 1

    print("nrois = %d" % (nrois))

    # === Load the matrix from the mat file === #
    check_file(inMatFN)
    condat = scipy.io.loadmat(inMatFN)
    
    assert(condat.keys().count("mn_cmat") == 1)
    assert(condat.keys().count("sprois") == 1)

    trois = deepcopy(condat["sprois"])
    trois = trois[0]
    assert(len(trois) == nrois)

    for (i0, t_roi) in enumerate(trois):
        t_str_roi = str(trois[i0])
        trois[i0] = t_str_roi.replace("[u'", "").replace("']", "")\
                    .replace("lh_", "L ").replace("rh_", "R ")
    trois = list(trois)

    idxr = []
    for (i0, t_roi) in enumerate(rois):
        idxr.append(trois.index(t_roi))
    trois = np.array(trois)

    tcon = deepcopy(condat["mn_cmat"][grp])
    mn_con = tcon[0][0]
#    mn_con = np.mean(tcon, axis=2)

    mn_con = mn_con[idxr, :]
    mn_con = mn_con[:, idxr]


    # == Set the self-connetions to zero == #
    for i0 in range(nrois):
        mn_con[i0][i0] = 0.0

    # === === #
    node_order = list(rois)
    node_angles = circular_layout(rois, node_order, start_pos=0)

    
    if np.isnan(vmax):
        vmax = np.max(mn_con)
        print("vmax = %.1f" % vmax)
    # con = np.random.rand(nrois, nrois) # DEBUG
    
    plot_connectivity_circle(mn_con, rois, node_angles=node_angles, 
                             facecolor="w", textcolor="k",
                             node_colors=roi_clrs,
                             colormap="binary",
                             vmax=vmax,
                             fontsize=12,
                             title="Connectivity matrix: %s - %s" % (grp, hemi))

    # === Save to tif file === #
    figFN = os.path.join(FIG_DIR, "conn_mat_circle_%s_%s.png" % (grp, hemi))
    pl.savefig(figFN, faceColor="w", format="png", dpi=200)
    check_file(figFN)
    print("INFO: Saved to image file: %s" % (figFN))

    pl.show()
