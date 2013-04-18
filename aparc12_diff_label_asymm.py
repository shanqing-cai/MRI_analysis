#!/usr/bin/python

import os
import sys
import argparse
import numpy as np
from subprocess import Popen, PIPE
from scipy import stats

from get_qdec_info import get_qdec_info
from scai_utils import *

LABEL_BASE_DIR = "/users/cais/STUT/analysis/aparc12_tracts_2"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Calculate the volume asymmetry index (a la Foundas et al.) of a certain aparc12 cortical label, in the diffusion space")
    ap.add_argument("subjsListFN", help="Subject list file name (e.g., /users/cais/STUT/FSDATA/qdec/subjsList_DTI.txt")
    ap.add_argument("label", help="aparc12 label name (e.g., PT)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    subjsListFN = args.subjsListFN
    label = args.label

    check_file(subjsListFN)

    f_subjsList = open(subjsListFN, 'r')
    t_subjs = f_subjsList.read().split('\n')
    f_subjsList.close()

    t_subjs = remove_empty_strings(t_subjs)

    sIDs = {"PFS": [], "PWS": []}
    for (i0, t_sID) in enumerate(t_subjs):
        t_grp = get_qdec_info(t_sID, "diagnosis")
        sIDs[t_grp].append(t_sID)
    
    grps = sIDs.keys()
    hemis = ["lh", "rh"]

    asymmIdx = {}
    for (i0, grp) in enumerate(grps):
        asymmIdx[grp] = [np.nan] * len(sIDs[grp])
        for (i1, sID) in enumerate(sIDs[grp]):
            sDir = os.path.join(LABEL_BASE_DIR, sID)
            check_dir(sDir)

            vols = [0] * len(hemis)
            for (j0, hemi) in enumerate(hemis):
                diff_mask_fn = os.path.join(sDir, \
                               "aparc12_%s_%s.diff.nii.gz" % (hemi, label))
                check_file(diff_mask_fn)

                [sout, serr] = Popen(["fslstats", diff_mask_fn, "-V"], 
                                     stdout=PIPE, stderr=PIPE).communicate()
                if len(serr) > 0:
                    raise Exception, "ERROR occurred during fslstats on subject %s, hemisphere %s: %s" % (sID, hemi, serr)

                vols[j0] = float(sout.split(' ')[0])          

            if vols.count(0) > 1:
                raise Exception, "ERROR occurred during loading data from subject %s" % sID
          
            iLeft = hemis.index("lh")
            iRight = hemis.index("rh")
            asymmIdx[grp][i1] = (vols[iLeft] - vols[iRight]) / \
                                (vols[iLeft] + vols[iRight]) * 2.0

        if asymmIdx[grp].count(np.nan) > 0:
            raise Exception, "Unable to calculate asymmIdx for some of the subjects in group %s" % grp
        
        print("Group %s (N=%d): mean = %f; ste = %f" % \
              (grp, len(asymmIdx[grp]), \
               np.mean(asymmIdx[grp]), \
               np.std(asymmIdx[grp]) / np.sqrt(len(asymmIdx[grp]))))
                

    (ttest_t, ttest_p) = stats.ttest_ind(asymmIdx["PWS"], asymmIdx["PFS"])
    print("t-test: t = %f; p = %f" % (ttest_t, ttest_p))
