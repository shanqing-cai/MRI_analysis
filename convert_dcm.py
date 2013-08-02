#!/speechlab/software/EPD7/epd-7.1-1-rh5-x86_64/bin/python

import os
import sys
import glob
import argparse
from scai_utils import *

help_doc = "Convert raw DICOM images to nii format. (scai)"

# === Heuristics === #
heur = [["functionalsparse", "bold"], 
        ["DIFFUSIONHIGHRES10Min", "diffusion"], 
        ["AAHeadScout32", "aascout"],
        ["tflmghmultiecho1mmiso", "T1"],
        ["localizer", "localizer"]]
        

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=help_doc)
    ap.add_argument("inputDir", type=str, \
                    help="Input raw DICOM data directory")
    ap.add_argument("outputDir", type=str, \
                    help="Output directory (for saving the .nii.gz files)")
    ap.add_argument("sID", type=str, \
                    help="Subject ID (e.g., ANS_M01)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    check_dir(args.inputDir)
    check_dir(args.outputDir, bCreate=True)

    # Create temporary directory #
    tmpDir = os.path.join(args.outputDir, "tmp_convert_dcm")
    check_dir(args.outputDir, bCreate=True)    
    delete_file_if_exists(tmpDir, recursive=True)
    check_dir(tmpDir, bCreate=True)
    

    
    dcm2nii_cmd = "dcm2nii -a n -o %s %s" \
                  % (os.path.abspath(tmpDir), \
                     os.path.join(os.path.abspath(args.inputDir), "*"))

    saydo(dcm2nii_cmd)
    sys.exit(0)

    # Apply heuristics to move and rename files #
    niis = glob.glob(os.path.join(tmpDir, "*.nii.gz"))
    niis.sort()
    for (i0, nii) in enumerate(niis):
        [nfp, nfn] = os.path.split(nii)

        imgStr = ""
        imgType = ""
        for (i1, t_item) in enumerate(heur):
            if (t_item[0] in nfn):
                imgStr = t_item[0]
                imgType = t_item[1]
                break

        if imgType == "":
            raise Exception, \
                  "Unrecognized image type for image file name %s" % nii
        
        # == Parse file name == #
        assert(nfn.endswith(".nii.gz"))
        assert(nfn.count("_") == 1)
        ius = nfn.index("_")
        assert(ius >= 8)
        timeStr = nfn[ius - 8 : ius + 7]
        assert(nfn.count(timeStr) == 1)

        prefix = nfn[:ius - 8]
        nfn1 = nfn.replace(timeStr, "")
        nfn2 = nfn1[len(prefix) :]
        suffix = nfn2.replace(imgStr, "").replace(".nii.gz", "")

        subDir = os.path.join(args.outputDir, imgType)
        check_dir(subDir, bCreate=True)

        if prefix == "":
            imgFN = os.path.join(subDir, \
                                  "%s_%s_%s.nii.gz" \
                                  % (args.sID, imgType, suffix))
        else:
            imgFN = os.path.join(subDir, \
                                  "%s_%s_%s_%s.nii.gz" \
                                  % (args.sID, imgType, suffix, prefix))

        saydo("cp %s %s" % (nii, imgFN))
        check_file(imgFN)


    # === Remove temporary directory === #
    saydo("rm -rf %s" % tmpDir)
