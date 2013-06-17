#!/usr/bin/python

import os
import sys
import argparse
import glob
import numpy as np

from scai_utils import *
from get_qdec_info import get_qdec_info

tractSegDir = "/users/cais/STUT/analysis/tractseg_aparc12/"
TRACULA_DIR = "/users/cais/STUT/analysis/dti2/tracula"
FNIRT_DIR = '/users/cais/STUT/analysis/nipype/T1_fnirt'
MNI152_TEMPLATE_FN = '/usr/share/fsl/data/standard/MNI152_T1_2mm.nii.gz'
L2_DIR = "/users/cais/STUT/analysis/tract_subcort_conn"
DESIGN_DIR = "/users/cais/STUT/analysis/design"

MATLAB_BIN = '/software/matlab2009a/bin/matlab'
con_01 = '/users/cais/STUT/analysis/design/con_01'


regions = ['Prefrontal', 'Premotor', 'Precentral', 'Postcentral', \
           'PPC', 'Occipital', 'Temporal']

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Group-level analysis of subcortical-cortical DTT connectivity")
    ap.add_argument("scSeed", help="Subcortical seed")
    ap.add_argument("--ccstop", dest="bCCStop", \
                    action="store_true", help="Use corpus callosum stop mask")
    ap.add_argument("--rerun", dest="bRerun", \
                    action="store_true", help="Rerun time-consuming steps")
    
    # === Process input arguments === #
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    
    scSeed = args.scSeed
    bCCStop = args.bCCStop
    bRerun = args.bRerun

    assert(scSeed.count(".") == 1)
    hemi = scSeed.split(".")[0]
    cSeed = scSeed.split(".")[1]


    # === Find subject IDs and group labels === #
    check_dir(tractSegDir)

    ds = glob.glob(os.path.join(tractSegDir, "*"))
    ds.sort()

    sIDs = []
    bPWS = []

    for (i0, t_fn) in enumerate(ds):
        (t_path, t_sID) = os.path.split(t_fn)
        sIDs.append(t_sID)
        
        if get_qdec_info(t_sID, "diagnosis") == "PWS":
            bPWS.append(1)
        else:
            bPWS.append(0)

    sIDs = np.array(sIDs)
    bPWS = np.array(bPWS)
    
    # === Locate the d2a files of each subject  === #
    d2as = []
    FAimgs = []
    warps = []
    for (i0, t_sID) in enumerate(sIDs):
        t_d2a = os.path.join(TRACULA_DIR, t_sID, "dmri", "xfms", "d2a.mat")
        check_file(t_d2a)
        d2as.append(t_d2a)

        t_FAimg = os.path.join(TRACULA_DIR, t_sID, "dmri", "dtifit_FA.nii.gz")
        check_file(t_FAimg)
        FAimgs.append(t_FAimg)

        s_fnirt_warp = os.path.join(FNIRT_DIR, t_sID, "T1_warp.nii.gz")
        check_file(s_fnirt_warp)
        warps.append(s_fnirt_warp)

    check_file(MNI152_TEMPLATE_FN)

    # === Transform the regMean files to subjects' structual spaces === #
    if bCCStop:
        scDirName = scSeed + "_ccStop"
    else:
        scDirName = scSeed
        
    regMean_ts = []
    for i0 in range(len(regions)):
        regMean_ts.append([])

    for (i0, t_sID) in enumerate(sIDs):
        sDir = os.path.join(tractSegDir, t_sID, scDirName)
        check_dir(sDir)

        for (i1, t_region) in enumerate(regions):
            regMean_d = os.path.join(sDir, \
                                     "tract_regionMean_%d.nii.gz" % (i1))
            check_file(regMean_d)
            
            regMean_t = os.path.join(sDir, \
                                     "tract_regionMean_%d.MNI152.nii.gz" % (i1))

            if bRerun or not os.path.isfile(regMean_t):
                warpCmd = 'applywarp --ref=' + MNI152_TEMPLATE_FN \
                          + ' --in=' + regMean_d \
                          + ' --warp=' + warps[i0] \
                          + ' --premat=' + d2as[i0] \
                          + ' --out=' + regMean_t
                saydo(warpCmd)

            check_file(regMean_t)
            regMean_ts[i1].append(regMean_t)

    # === Concatenate subject files into 4D files === 
    # === in prep for group-level analysis === #
    check_dir(L2_DIR)
    regMean_4ds = []
    L2SeedDir = os.path.join(L2_DIR, scDirName)
    check_dir(L2SeedDir, bCreate=True)

    for (i0, t_region) in enumerate(regions):
        regMean_4ds.append(os.path.join(L2SeedDir, "merged_%d.nii.gz" % i0))
        mergeCmd= "fslmerge -t %s " % regMean_4ds[-1]

        for (i1, t_fn) in enumerate(regMean_ts[i0]):
            mergeCmd += "%s " % t_fn
        
        if bRerun or not os.path.isfile(regMean_4ds[-1]):
            saydo(mergeCmd)
            check_file(regMean_4ds[-1])

    # === Generate the design matrix  === #
    X_line = 'X = ['
    
    for (i0, s) in enumerate(sIDs):
        if bPWS[i0] == 1:
            X_line += '1; '
        else:
            X_line += '-1; '

    X_line = X_line[:-2]
    X_line += '];\nX = [ones(' + str(len(sIDs)) + ', 1), X];'

    check_dir(DESIGN_DIR)
    X_mat = os.path.join(DESIGN_DIR, "X_subcort_conn.mat")
    X_line += '\nsave(\'' + X_mat + '\', \'X\', \'-v4\')\n'

    mScriptGenX = os.path.join(DESIGN_DIR, 'genX_subcort_conn.m')
    genX_f = open(mScriptGenX, 'w')
    genX_f.write(X_line)
    genX_f.close()
    check_file(mScriptGenX)

    matlabCmd = MATLAB_BIN + ' -nosplash -nodesktop -r \'run ' \
                + mScriptGenX + '; exit\''
    saydo(matlabCmd)
    check_file(X_mat)

    # === Perform between-group comparisons === #
    check_file(con_01)
    
    for (i0, t_region) in enumerate(regions):
        t_con_01_dir = os.path.join(L2SeedDir, "bgc_%d" % i0)
        fitCmd = 'mri_glmfit --y ' + regMean_4ds[i0] \
                 + ' --X ' + X_mat \
                 + ' --C ' + con_01 \
                 + ' --glmdir ' + t_con_01_dir

        saydo(fitCmd)
        check_dir(t_con_01_dir)

        sig_fn = os.path.join(t_con_01_dir, "con_01", "sig.mgh")
        check_file(sig_fn)
    
        
