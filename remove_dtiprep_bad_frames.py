#!/usr/bin/python

import os
import sys
import glob
import argparse
import tempfile
from scai_utils import *
from dtiprep_qc_xml import get_bad_frames

DATA_DIR = "/users/cais/STUT/DATA"
DTIPREP_DIR = "/users/cais/STUT/analysis/dti"

DATA_DIR_RHY = "/users/cais/RHY/DATA"
DTIPREP_DIR_RHY = "/users/cais/RHY/analysis/dwi"

DEFAULT_N_FRAMES = 70

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Remove frames that are determined as bad by DTIPrep QC from the original nifti DWI")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("--nFrames", type=int, dest="nFrames", \
                    default=DEFAULT_N_FRAMES, \
                    help="Number of volumes (default = %d)" % DEFAULT_N_FRAMES)
    ap.add_argument("--RHY", dest="bRHY", action="store_true", \
                    help="RHY project, instead of the default STUT project")

    if len(sys.argv) == 1:
        ap.print_help()
        rsys.exit(1)

    args = ap.parse_args()
    sID = args.sID

    if args.bRHY:
        DATA_DIR = DATA_DIR_RHY
        DTIPREP_DIR = DTIPREP_DIR_RHY
    
    # Locate the original diffusion data
    sDataDir = os.path.join(DATA_DIR, sID)
    check_dir(sDataDir)
    sOrigDiffDir = os.path.join(sDataDir, "diffusion")
    check_dir(sOrigDiffDir)

    if not args.bRHY:
        d0 = glob.glob(os.path.join(sOrigDiffDir, "%s_run??_??.nii.gz" % sID))
    else:
        d0 = glob.glob(os.path.join(sOrigDiffDir, \
                                    "%s_diffusion_*.nii.gz" % sID))

    d0.sort()
    if len(d0) == 0:
        raise Exception, "Failed to find diffusion series in directory %s" \
                         % d0
    
    if not args.bRHY:
        baseName = d0[-1].replace(".nii.gz", "")
    else:
        baseName = d0[0].replace(".nii.gz", "")
    
    dwi4dName = baseName + ".nii.gz"

    if not args.bRHY:
        bvecName = baseName + ".mghdti.bvecs"
        bvalName = baseName + ".mghdti.bvals"
    else:
        bvecName = baseName + ".bvec"
        bvalName = baseName + ".bval"

    check_file(bvecName)
    check_file(bvalName)

    print("dwi4dName = %s" % dwi4dName)
    print("bvalName = %s" % bvalName)
    print("bvecName = %s" % bvecName)

    # Locate the QC result XML
    sQCDir = os.path.join(DTIPREP_DIR, sID)
    check_dir(sQCDir)
    qcXML = os.path.join(sQCDir, "%s_XMLQCResult.xml" % sID)
    check_file(qcXML)
    
    badFrmIdx = get_bad_frames(qcXML, nFrames=args.nFrames)

    print("\nThere are %d bad frame(s) according to %s:" \
          % (len(badFrmIdx), qcXML))
    if len(badFrmIdx) > 0:
        print(badFrmIdx)


    # Eliminate the bad frames from the dwi4d
    qcedDir = os.path.join(sOrigDiffDir, "qced")
    check_dir(qcedDir, bCreate=True)
    (tpath, tbase) = os.path.split(baseName)
    qced4d = os.path.join(qcedDir, tbase + "_qced.nii.gz")

    if not args.bRHY:
        (foo, baseName0)= os.path.split(bvalName)
        qced_bvalName = os.path.join(qcedDir, \
                                 baseName0.replace(".mghdti", "_qced.mghdti"))
        (foo, baseName0) = os.path.split(bvecName)
        qced_bvecName = os.path.join(qcedDir, \
                                 baseName0.replace(".mghdti", "_qced.mghdti"))
    else:
        (foo, baseName0)= os.path.split(bvalName)
        qced_bvalName = os.path.join(qcedDir, \
                                 baseName0.replace(".bval", ".qced.bval"))
        (foo, baseName0) = os.path.split(bvecName)
        qced_bvecName = os.path.join(qcedDir, \
                                 baseName0.replace(".bvec", ".qced.bvec"))
    
    
    tmpName = tempfile.mktemp()
    (foo, tmpName) = os.path.split(tmpName)
    tmpDir = os.path.join(sOrigDiffDir, tmpName)
    check_dir(tmpDir, bCreate=True)

    if len(badFrmIdx) == 0:
        saydo("cp %s %s" % (dwi4dName, qced4d))
        check_file(qced4d)

        saydo("cp %s %s" % (bvalName, qced_bvalName))
        check_file(qced_bvalName)

        saydo("cp %s %s" % (bvecName, qced_bvecName))
        check_file(qced_bvecName)
    else:
        splitOutBase = os.path.join(tmpDir, 'vol')
        split_cmd = 'fslsplit %s %s -t' % (dwi4dName, splitOutBase)
        saydo(split_cmd)
        
        for t_bfi in badFrmIdx:
            badFrmFN = splitOutBase + "%.4d" % t_bfi + ".nii.gz"
            check_file(badFrmFN)
            saydo("rm -f %s" % badFrmFN)
            
        # Re-merge
        dvols = glob.glob(splitOutBase + "????.nii.gz")
        dvols.sort()
        if len(dvols) != args.nFrames - len(badFrmIdx):
            raise Exception, "Unexpected error: wrong number of frames after removal of QC-failure frames"

        merge_cmd = "fslmerge -t %s " % qced4d
        for t_vol in dvols:
            merge_cmd += "%s " % t_vol

        saydo(merge_cmd)
        check_file(qced4d)

        saydo("rm -r %s" % tmpDir)

        # Copy and fix the bvals
        bvals_f0 = open(bvalName, 'r')
        bvals_txt0 = bvals_f0.read().split('\n')
        bvals_txt0 = remove_empty_strings(bvals_txt0)
        bvals_f0.close()
        if len(bvals_txt0) != args.nFrames:
            raise Exception, \
                  "Erroneous number of lines in bvals file: %s" % bvals_txt0
        
        bvals_f1 = open(qced_bvalName, "w")
        for (i0, t_line) in enumerate(bvals_txt0):
            if badFrmIdx.count(i0) == 0:
                bvals_f1.write("%s\n" % t_line)
        bvals_f1.close()
        check_file(qced_bvalName)

        # Copy and fix the bvecs
        bvecs_f0 = open(bvecName, 'r')
        bvecs_txt0 = bvecs_f0.read().split('\n')
        bvecs_txt0 = remove_empty_strings(bvecs_txt0)
        bvecs_f0.close()
        if len(bvecs_txt0) != args.nFrames:
            raise Exception, \
                  "Erroneous number of lines in bvecs file: %s" % bvecs_txt0
        
        bvecs_f1 = open(qced_bvecName, "w")
        for (i0, t_line) in enumerate(bvecs_txt0):
            if badFrmIdx.count(i0) == 0:
                bvecs_f1.write("%s\n" % t_line)
        bvecs_f1.close()
        check_file(qced_bvecName)
        

