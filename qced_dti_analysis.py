#!/usr/bin/python

import os
import sys
import glob
import argparse
from scai_utils import *

REQ_FSL_VER = "5.0"

DATA_DIR = "/users/cais/STUT/DATA"
DTIPREP_DIR = "/users/cais/STUT/analysis/dti"
DTI_BASE = "/users/cais/STUT/analysis/dti2"
EDDY_CORRECT_SPLINE_BIN = "/users/cais/STUT/scripts/eddy_correct_spline"
#FDT_ROTATE_BVECS_BIN = "/users/cais/STUT/scripts/fdt_rotate_bvecs"
FDT_ROTATE_BVECS_BIN = "/users/cais/STUT/scripts/rotate_bvecs.py"
TRACULA2_BIN = "/users/cais/STUT/scripts/tracula2.py"
BEDP_SCRIPT = "/users/cais/STUT/scripts/trac-all_bedp_3.sh"

DATA_DIR_RHY = "/users/cais/RHY/DATA"
DTIPREP_DIR_RHY = "/users/cais/RHY/analysis/dwi"
DTI_BASE_RHY = "/users/cais/RHY/analysis/dwi"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run DWI analysis on qced 4d data (see remove_dtiprep_bad_frames.py")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("--rerun", dest="bForceRerun", action="store_true", \
                    help="Force rerun finished steps")
    ap.add_argument("--RHY", dest="bRHY", action="store_true", \
                    help="Project RHY, instead of the default project STUT")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    sID = args.sID
    bForceRerun = args.bForceRerun

    if args.bRHY:
        DATA_DIR = DATA_DIR_RHY
        DTIPREP_DIR = DTIPREP_DIR_RHY
        DTI_BASE = DTI_BASE_RHY

    # Check the version of FSL
    env = os.environ
    if env["FSLDIR"].count(REQ_FSL_VER) == 0:
        raise Exception, "It appears that a version of FSL other than %s is being used" % REQ_FSL_VER
    
    # Locate the original 4d input and accompanying bvals and bvecs files
    qcedDir = os.path.join(DATA_DIR, sID, "diffusion", "qced");
    check_dir(qcedDir)

    if not args.bRHY:
        d0 = glob.glob(os.path.join(qcedDir, "%s_run??_??_qced.nii.gz" % sID))
    else:
        d0 = glob.glob(os.path.join(qcedDir, \
                                    "%s_diffusion_*_qced.nii.gz" % sID))

    if len(d0) != 1:
        raise Exception, "Not exactly one 4D series found in directory: %s" % \
                         qcedDir

    dwi4d = d0[0];

    if not args.bRHY:
        bvalsFN = dwi4d.replace(".nii.gz", ".mghdti.bvals")
        bvecsFN = dwi4d.replace(".nii.gz", ".mghdti.bvecs")
    else:
        bvalsFN = dwi4d.replace("_qced.nii.gz", ".qced.bval")
        bvecsFN = dwi4d.replace("_qced.nii.gz", ".qced.bvec")

    check_file(bvalsFN)
    check_file(bvecsFN)
    
    # Execute eddy_correct_spline
    outBaseDir = os.path.join(DTI_BASE, sID)
    check_dir(outBaseDir, bCreate=True)
    
    ec4d = os.path.join(outBaseDir, "ecdwi.nii.gz")
    ec_cmd = "%s %s %s 0" % (EDDY_CORRECT_SPLINE_BIN, \
                           dwi4d, ec4d)
    ecclog = os.path.join(outBaseDir, "ecdwi.ecclog")


    b_ecDone = os.path.isfile(ec4d) and os.path.isfile(ecclog)
    if not b_ecDone or bForceRerun:
        saydo("rm -f %s" % ecclog)
        saydo(ec_cmd)
        check_file(ec4d)
        check_file(ecclog)
    else:
        if not bForceRerun:
            print("Skipping step eddy_correct_splint (already done)")

    # Execute fdt_rotate_bvecs
    rotBvecsFN = os.path.join(outBaseDir, "rotated.bvecs")
    rotBvalsFN = os.path.join(outBaseDir, "rotated.bvals")
    rot_cmd = "%s %s %s %s -v" % (FDT_ROTATE_BVECS_BIN, \
                               bvecsFN, rotBvecsFN, ecclog)

    b_rotDone = os.path.isfile(rotBvecsFN) and os.path.isfile(rotBvalsFN)
    if not b_rotDone or bForceRerun:        
        saydo(rot_cmd)
        check_file(rotBvecsFN)

        saydo("cp %s %s" % (bvalsFN, rotBvalsFN))
        check_file(rotBvalsFN)
    else:
        if not bForceRerun:
            print("Skipping step fdt_rotate_bvecs (already done)")
            
    # Execute tracula2.py
    tracula2_cmd = "%s %s prep" % (TRACULA2_BIN, sID)
    if args.bRHY:
        tracula2_cmd += " --RHY"
    
    # Determine if tracula prep has already finished
    tracBase = os.path.join(DTI_BASE, "tracula")
    sTracDir = os.path.join(tracBase, sID)
    bTraculaPrepDone = \
        os.path.isfile(os.path.join(sTracDir, "dmri", "bvals")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "bvecs")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_FA.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_L1.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_L2.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_L3.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_MD.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_MO.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_S0.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_V1.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_V2.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dtifit_V3.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "dwi.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "lowb.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "lowb_brain.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "brain_anat_mni.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dmri", "brain_anat.nii.gz")) and \
        os.path.isfile(os.path.join(sTracDir, "dlabel", "diff", "lowb_brain_mask.nii.gz"))
        

    if not bTraculaPrepDone or bForceRerun:
        saydo(tracula2_cmd)
    else:
        print("Skipping step tracula2.py prep (already done)")

    # Run bedpostx (trac-all_bedp_3.sh)
    check_file(BEDP_SCRIPT)
    bedp_cmd = "%s %s" % (BEDP_SCRIPT, sID)
    
    # Determine if bedp is already done
#    bvecs1 = 
    
        
    
    
