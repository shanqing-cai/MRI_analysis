#!/usr/bin/python

import os
import sys
import tempfile
from subprocess import Popen, PIPE

from scai_utils import *
from aparc12 import get_aparc12_cort_rois


FSAV_DIR = "/users/cais/STUT/FSDATA/fsaverage2"
APARC12_VOL = "/users/cais/STUT/FSDATA/fsaverage2/mri/aparcJT+aseg.mgz"
APARC12_CTAB = "/users/cais/STUT/slaparc_550.ctab"

hemis = ["lh", "rh"]

if __name__ == "__main__":
    check_dir(FSAV_DIR)
    check_file(APARC12_VOL)
    check_file(APARC12_CTAB)

    # Read color table
    (roi_nums, roi_names) = read_ctab(APARC12_CTAB)

    # Get the set of ROIs to process
    proc_rois_0 = get_aparc12_cort_rois(bSpeech=False)

    # Duplex the rois in to both hemispheres and get their numbers
    proc_roi_names = []
    proc_roi_nums = []

    for (i0, themi) in enumerate(hemis):
        for (i1, troi) in enumerate(proc_rois_0):
            troi_name = themi + "_" + troi
            proc_roi_names.append(troi_name)
        
            if roi_names.count(troi_name) == 0:
                raise Exception, "Cannot find ROI %s in color table %s" \
                                 % (troi_name, APARC12_CTAB)
            
            proc_roi_nums.append(roi_nums[roi_names.index(troi_name)])
    
    assert len(proc_roi_names) == len(proc_roi_nums)

    # Open file for writing 
    mriDir = os.path.join(FSAV_DIR, "mri")
    check_dir(mriDir)
    resFN = os.path.join(mriDir, "aparc12_roi_coords.txt")
    resF = open(resFN, "wt")

    # Iterate through all ROIs
    for (i0, t_roi_name) in enumerate(proc_roi_names):
        binOut = tempfile.mktemp() + ".nii.gz"
        t_roi_num = proc_roi_nums[i0]
        
        binCmd = "mri_binarize --i %s --match %d --o %s" \
                 % (APARC12_VOL, t_roi_num, binOut)
        saydo(binCmd)
        check_file(binOut)
        
        (sout, serr) = Popen(["fslstats", binOut, "-C"], \
                             stdout=PIPE, stderr=PIPE).communicate()
        saydo("rm -f %s" % binOut)

        if len(serr) > 0:
            raise Exception, "Error occurred during fslstats on the binarized image for ROI %s (roi_code=%d)" % (t_roi_name, t_roi_num)
        
        t_cog = [-1] * 3
        for i1 in range(3):
            t_cog[i1] = float(sout.split(' ')[i1])
            
        assert(t_cog.count(-1) == 0)

        resF.write("%s %d %f %f %f\n" % \
                   (t_roi_name, t_roi_num, t_cog[0], t_cog[1], t_cog[2]))

    resF.close()
    check_file(resFN)
    
    print("Results (center-of-gravity) data saved to file: %s" % resFN)
