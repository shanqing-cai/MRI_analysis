#!/usr/bin/python

import os
import sys
import glob
import string
import argparse
import numpy as np
from subprocess import Popen, PIPE
from scipy.io import savemat
from scai_utils import *


## Config: paths
FNIRT_DIR = '/users/cais/STUT/analysis/nipype/normalize'
FSDATA_DIR = '/users/cais/STUT/FSDATA'
TRACULA_DIR = '/users/cais/STUT/tracula/'
TRACULA_DIR_DTIPREP = '/users/cais/STUT/analysis/dti2/tracula'
MNI152_TEMPLATE_FN = "/usr/share/fsl/data/standard/MNI152_T1_2mm.nii.gz"
GROUP_MNI_DIR = "/users/cais/STUT/analysis/funcloc_trac"
DATA_DIR ='/users/cais/STUT/DATA'
ROI_TRACK_BASE = '/users/cais/STUT/ROI_TRACTS/'
CTAB_FN = '/users/cais/STUT/slFRS17.ctab'

VALID_MEAS = ["FA", "L1", "RD"];

MULT_FACTOR = 1000.0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate whie-matter FA tables of different depths based on aparc12")
    ap.add_argument("sID", type=str, help="Subject ID")
    ap.add_argument("xmm", type=str, help="White-matter projection depth (e.g., 3)")
    ap.add_argument("--dtiprep", dest="bDTIPrep", action="store_true", \
                    help="Use the DTIPrep preprocessed data")
    ap.add_argument("--meas", dest="meas", default="FA", \
                    help="Specify diffusion-tensor measure type {FA, L1, RD}. RD is defined as the arithmetic mean of L2 and L3")
    ap.add_argument("--rerunbbr", dest="bRerunBBR", action="store_true", \
                    help="Do not use the exisiting bbr mat; but rerun BBR to generate a new mat")
    ap.add_argument("--rerun", dest="bRerun", action="store_true", \
                    help="Force rerun time-consuming steps")


    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    sID = args.sID
    xmm = int(args.xmm)
    bDTIPrep = args.bDTIPrep
    meas = args.meas
    bRerunBBR = args.bRerunBBR
    bRerun = args.bRerun

    # Input sanity check
    if VALID_MEAS.count(meas) == 0:
        raise Exception, "Unrecognized diffusion tensor measure: %s" % meas

    if bDTIPrep:
        TRACULA_DIR = TRACULA_DIR_DTIPREP

    trac_dir = os.path.join(TRACULA_DIR, sID)
    if os.path.isdir(trac_dir):
        print('TRACULA directory = ' + trac_dir)
    else:
        sys.exit('ERROR: DATA directory not found: %s'%trac_dir)

    
    if meas == 'FA':
        dtmeas_fn = os.path.join(trac_dir, 'dmri', 'dtifit_FA.nii.gz')
    else:
        dtmeas_fn = os.path.join(trac_dir, 'dmri', 'dtifit_' + meas + '.nii.gz')
    
    if not os.path.isfile(dtmeas_fn):
        raise Exception, "Cannot find diffusion-tensor measure (%s) file: %s" \
                         % (meas, dtmeas_fn)
    else:
        print("dtmeas_fn = %s"%dtmeas_fn)

    dtmeas_anat_fn = os.path.join(trac_dir, 'dmri', 'dtifit_' + meas + '.anat.nii.gz')
    brain_mgz = os.path.join(FSDATA_DIR, sID, 'mri', 'brain.mgz')
    brain_ngz = os.path.join(FSDATA_DIR, sID, 'mri', 'brain.nii.gz')
    if os.path.isfile(brain_ngz) and (not bRerun):
        print("INFO: brain.nii.gz already exists: %s"%brain_ngz)
    else:
        cvt_cmd = 'mri_convert ' + brain_mgz + ' ' + brain_ngz
        saydo(cvt_cmd)
    
    

    xfm_diff2anat_mat = os.path.join(trac_dir, 'dmri', 'xfms', 'diff2anatorig.bbr.mat')
    #xfm_diff2anat_mat = os.path.join(trac_dir, 'dmri', 'xfms', 'anatorig2diff.bbr.mat')
#    xfm_anat2diff_mat = os.path.join(trac_dir, 'dmri', 'xfms', 'anatorig2diff.bbr.mat')
    #if not os.path.isfile(xfm_anat2diff_mat):
    #    raise Exception, "Cannot find anatorig-to-diff transform mat: %s"%xfm_anat2diff_mat

    segfn = os.path.join(DATA_DIR, sID, "aparc12_wm%dmm.nii.gz" % (xmm))
    if not os.path.isfile(segfn):
        raise Exception, "Cannot find aparc12 wm seg file: %s. Use aparc12_wm.py to generate it first."%(segfn)
    print("INFO: segfn = %s" % segfn)

#    segdiff_fn = segfn.replace(".nii.gz", ".diff.nii.gz")
#    a2d_xfm_cmd = "flirt -in %s -ref %s -applyxfm -init %s -out %s" \
#                  % (segfn, dtmeas_fn, xfm_anat2diff_mat, segdiff_fn)
        #xfm_diff2anat_mat = os.path.join(trac_dir, 'dmri', 'xfms', 'new_diff2anatorig.bbr.mat')
        #xfm_diff2anat_mat = os.path.join(trac_dir, 'dmri', 'xfms', 'new_diff2anatorig.bbr.mat')
        #bbr_cmd = "bbregister --s %s --mov %s --init-spm --dti --reg %s --fslmat %s"\
        #          % (sID, )
        
    if not os.path.isfile(xfm_diff2anat_mat):
        sys.exit('ERROR: diff2anatorig.bbr.mat not found')
    else:
        print('diff2anatorig xfm_mat = ' + xfm_diff2anat_mat)


    if os.path.isfile(dtmeas_anat_fn) and (not bRerun):
        print("INFO: anatomical-space volume already exists: %s"%dtmeas_anat_fn)
    else:
        flirt_cmd = 'flirt -in ' + dtmeas_fn + ' -ref ' + brain_ngz + \
                    ' -applyxfm -init ' + xfm_diff2anat_mat + \
                    ' -out ' + dtmeas_anat_fn
        saydo(flirt_cmd)

    # Multiply by a factor for L1 and RD (for precision issues)
    if meas == "L1" or meas == "RD":
        bMult = True
        mult_cmd = "fslmaths %s -mul %f %s" \
                   % (dtmeas_anat_fn, MULT_FACTOR, dtmeas_anat_fn)
        saydo(mult_cmd)
        check_file(dtmeas_anat_fn)
    else:
        bMult = False

    # ------ Run mri_segstats to get stats info ------ #
    sumfn = os.path.join(trac_dir, 'dmri', \
                         "aparc12_2_%s_wm%dmm.sum"%(meas, xmm))
    if os.path.isfile(sumfn) and (not bRerun):
        print("INFO: sum file already exists: %s. Not rerunning (see --rerun option)."%sumfn)
    else:
        segfn = os.path.join(DATA_DIR, sID, "aparc12_wm%dmm.nii.gz"%(xmm))
        if not os.path.isfile(segfn):
            raise Exception, "Cannot find aparc12 wm seg file: %s. Use aparc12_wm.py to generate it first."%(segfn)
        else:
            print("segfn = %s"%(segfn))
        cmd = "mri_segstats --seg %s --in %s --sum %s"\
              %(segfn, dtmeas_anat_fn, sumfn)
        saydo(cmd)

    

    # ------ Parse the sum file ------ #    
    sumf = open(sumfn, 'r')
    stxt = sumf.read().split('\n')
    sumf.close()

    # Read color table
    if not os.path.isfile(CTAB_FN):
        raise Exception, "Cannot open color table (ctab): %s"%(CTAB_FN)
    (ct_nums, ct_names) = read_ctab(CTAB_FN)

    # Go to the title line
    bFoundTitleLine = False
    for (i0, t_line) in enumerate(stxt):
        if t_line.startswith('# ColHeaders'):
            bFoundTitleLine = True
            break
    if (not bFoundTitleLine):
        raise Exception, "Cannot parse sum file: %s"%(sumfn)
    stxt = stxt[i0 + 1:]

    rois = []
    meanval = []
    for (i0, t_line) in enumerate(stxt):
        if len(t_line) == 0:
            continue
        
        t_line = t_line.replace(' ', '\t')
        t_items = t_line.split('\t')
        while t_items.count('') > 0:
            t_items.remove('')
        
        t_rnum = int(t_items[1])
        if t_rnum < 3000 or t_rnum > 4999:
            continue
        print("ROI index = %d"%t_rnum)

        if t_rnum / 1000 == 3:
            hemi = "lh"
        elif t_rnum / 1000 == 4:
            hemi = "rh"
        else:
            raise Exception, "Unexpected ROI index: %d"%t_rnum

        t_rnum = t_rnum - (t_rnum / 1000) * 1000
        
        # Find out ROI name
        bFound = False
        for (i1, t_num) in enumerate(ct_nums):
            if t_num == t_rnum:
                bFound = True
                break
        if not bFound:
            print("WARNING: cannot find ROI number %d in color table."%t_rnum)
            continue
        t_roi = "%s_%s"%(hemi, ct_names[i1])
        print("ROI name = %s"%(t_roi))

        rois.append(t_roi)
        if not bMult:
            meanval.append(float(t_items[5]))
        else:
            meanval.append(float(t_items[5]) / MULT_FACTOR)

    
    # Save to mat file
    meanval = np.array(meanval)
    
    faData = {"rois": rois, "mean%s" % meas: meanval}
    if not bDTIPrep:
        matfn = os.path.join(DATA_DIR, sID, "aparc12_2_%s_wm%smm.mat"%(meas, xmm))
    else:
        matfn = os.path.join(DATA_DIR, sID, "aparc12_dtiprep2_%s_wm%smm.mat"%(meas, xmm))
    savemat(matfn, faData)
    
    if not os.path.isfile(matfn):
        raise Exception, "It appears that savemat failed to save FA data to .mat file: %s"%matfn
    else:
        print("Saved %s data to .mat file: %s"%(meas, matfn))
