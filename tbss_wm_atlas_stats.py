#!/usr/bin/python

import os
import sys
import glob
import argparse
import tempfile
import numpy as np
from scipy.io import *
from scipy import stats
from subprocess import Popen, PIPE

from scai_utils import *
from get_qdec_info import get_qdec_info
from read_xml_labels import read_xml_labels

atlas_label_fn = \
    "/usr/share/fsl/5.0/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz"
atlas_label_xml = \
    "/usr/share/fsl/5.0/data/atlases/JHU-labels.xml"

P_THRESH_UNC = 0.05

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Get stats (e.g., average FA) from in atlas-defined WM regions in TBSS-aligned diffusion-tensor images")
    ap.add_argument("tbssDir", help="Base TBSS directory (e.g., /users/cais/STUT/analysis/dt_tbss_dtiprep2)")

    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
        
    # === Parse input arguments === #
    args = ap.parse_args()
    tbssDir = args.tbssDir

    # === Input sanity check === #
    check_dir(tbssDir)

    statsDir = os.path.join(tbssDir, "stats")
    check_dir(statsDir)

    origDir = os.path.join(tbssDir, "origdata")
    check_dir(origDir)

    check_file(atlas_label_fn)

    # === Read JHU-ICBM labels === #
    check_file(atlas_label_xml)
    
    labs = read_xml_labels(atlas_label_xml)
    
    # === Locate the all_FA image === #
    allFA = os.path.join(statsDir, "all_FA.nii.gz")
    check_file(allFA)
    
    # === Find out the subject IDs and their groups === #
    origDir = os.path.join(tbssDir, "origdata")
    check_dir(origDir)

    ds = glob.glob(os.path.join(origDir, "S??.nii.gz"))
    ds.sort()

    sIDs = []
    idxPWS = []
    idxPFS = []

    for (i0, d) in enumerate(ds):
        [tpath, tfn] = os.path.split(d)
        sID = tfn.replace(".nii.gz", "")
        sIDs.append(sID)
        
        if get_qdec_info(sID, "diagnosis") == "PWS":
            idxPWS.append(i0)
        elif get_qdec_info(sID, "diagnosis") == "PFS":
            idxPFS.append(i0)
        else:
            raise Exception, "Unrecognized diagnosis for subject %s: %s" % \
                             (sID, get_qdec_info(sID, "diagnosis"))

    # === Split the all_FA image, for later use by fslstats  === #
    splitBase = tempfile.mktemp()
    split_cmd = "fslsplit %s %s -t" % (allFA, splitBase)
    saydo(split_cmd)
    
    splitFNs = glob.glob(splitBase + "*.nii.gz")
    splitFNs.sort()

    if len(splitFNs) != len(sIDs):
        raise Exception, "Number of volumes in 4D series %s (%d) does not match the number of subjects in origdata (%d)" % \
                         (allFA, len(splitFNs), len(sIDs))

    # === Iterate through the WM labels and get the stats info === #
    labRes = {"labels": [], "meanFA": [], "tt_t": [], "tt_p": []}

    for (i0, lab) in enumerate(labs['name']):
        ind = labs['ind'][i0]
        
        if ind == 0:
            continue

        print("\nProcessing label #%d: %s\n" % (i0, lab))
        
        labRes["labels"].append(lab)
        labRes["meanFA"].append({"PWS": [], "PFS": []})

        tmpfn = tempfile.mktemp() + ".nii.gz"

        # == Binarize, get label mask == #
        bin_cmd = "mri_binarize --i %s --match %d --o %s" % \
                  (atlas_label_fn, ind, tmpfn)
        saydo(bin_cmd)
        check_file(tmpfn)

        # == Use fslstats to get the masked mean == #
        t_vals = [-1] * len(sIDs)
        
        for (i1, splitFN) in enumerate(splitFNs):            
            (sout, serr) = Popen(["fslstats", splitFN, "-k", tmpfn, "-m"], \
                                 stdout=PIPE, stderr=PIPE).communicate()

            if len(serr) > 0:
                raise Exception, \
                      "ERROR occurred during fslstats on %s" % splitFN

            t_vals[i1] = float(sout.split(' ')[0])
        
        t_vals = np.array(t_vals)
        labRes["meanFA"][-1]["PWS"] = t_vals[idxPWS]
        labRes["meanFA"][-1]["PFS"] = t_vals[idxPFS]

        (t, p) = stats.ttest_ind(labRes["meanFA"][-1]["PWS"], \
                                 labRes["meanFA"][-1]["PFS"])

        labRes["tt_t"].append(t)
        labRes["tt_p"].append(p)
        
        os.system("rm -f %s" % tmpfn)
    
    os.system("rm -f %s*" % splitBase)
    
    # === Save results to mat file === #
    resMatFN = "/users/cais/STUT/scripts/tbss_wm_atlas_stats.mat"
    os.system("rm -f %s" % resMatFN)
    savemat(resMatFN, labRes)
    check_file(resMatFN)

    print("\nINFO: Results saved to .mat file: %s" % resMatFN)

    # === Print results === #
    print("=== Significant results at P_THRESH_UNC = %f ===" % P_THRESH_UNC)
    for (i0, labName) in enumerate(labRes["labels"]):
        if labRes["tt_p"][i0] < P_THRESH_UNC:
            mean_PFS = np.mean(labRes["meanFA"][i0]["PFS"])
            mean_PWS = np.mean(labRes["meanFA"][i0]["PWS"])
            ste_PFS = np.std(labRes["meanFA"][i0]["PFS"]) / \
                      np.sqrt(len(idxPFS))
            ste_PWS = np.std(labRes["meanFA"][i0]["PWS"]) / \
                      np.sqrt(len(idxPWS))
            
            print("WM label [%s]:" % labName)
            print("\tPFS: mean = %f; SE = %f" % (mean_PFS, ste_PFS))
            print("\tPWS: mean = %f; SE = %f" % (mean_PWS, ste_PWS))
            print("\tt = %f; p = %f" % \
                  (labRes["tt_t"][i0], labRes["tt_p"][i0]))
