#!/usr/bin/python

import os
import sys
import argparse
import glob
import numpy as np

from get_qdec_info import get_qdec_info
from scai_utils import *

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="funcloc data subject demographic summary")
    ap.add_argument("inDir", help="funcloc base directory")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    inDir = args.inDir

    check_dir(inDir)
    
    d0 = glob.glob(os.path.join(inDir, "S??"))

    sIDs = []

    for (i0, t_d) in enumerate(d0):
        sIDs.append(os.path.split(t_d)[1])
        

    sIDs.sort()

    grps = []
    ages = []
    genders = []
    for (i0, t_sID) in enumerate(sIDs):
        if get_qdec_info(t_sID, "diagnosis") == "PWS":
            grps.append(1)
        else:
            grps.append(0)

        if get_qdec_info(t_sID, "gender") == "Male":
            genders.append(1)
        else:
            genders.append(0)

        ages.append(float(get_qdec_info(t_sID, "Age")))

    grps = np.array(grps)
    ages = np.array(ages)
    genders = np.array(genders)

    # Print summary
    ages_PWS = ages[np.nonzero(grps)[0]]
    print("PWS: N=%d, %dF%dM; Age min=%.1f, max=%.1f, mean=%.1f, SD=%.1f" % \
          (len(np.nonzero(grps)[0]), \
           len(np.intersect1d(np.nonzero(grps)[0], np.nonzero(1 - genders)[0])), \
           len(np.intersect1d(np.nonzero(grps)[0], np.nonzero(genders)[0])), \
           np.min(ages_PWS), np.max(ages_PWS), np.mean(ages_PWS), np.std(ages_PWS)))
    
    ages_PFS = ages[np.nonzero(1 - grps)[0]]
    print("Controls (PFS): N=%d, %dF%dM; Age min=%.1f, max=%.1f, mean=%.1f, SD=%.1f" % \
          (len(np.nonzero(1 - grps)[0]), \
           len(np.intersect1d(np.nonzero(1 - grps)[0], np.nonzero(1 - genders)[0])), \
           len(np.intersect1d(np.nonzero(1 - grps)[0], np.nonzero(genders)[0])), \
           np.min(ages_PFS), np.max(ages_PFS), np.mean(ages_PFS), np.std(ages_PFS)))
