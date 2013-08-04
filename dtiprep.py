#!/usr/bin/python

import os
import sys
import argparse
import glob
from scai_utils import *

DATA_DIR = "/users/cais/STUT/DATA"
DATA_DIR_RHY = "/users/cais/RHY/DATA"
FSDATA_DIR = "/users/cais/STUT/FSDATA"
# DTIPREP_PATH = "/software/DTIPrep/1.1.6"
DTIPREP_PATH = "/software/DTIPrep/130630"
OUTPUT_DIR = "/users/cais/STUT/analysis/dti"
OUTPUT_DIR_RHY = "/users/cais/RHY/analysis/dwi"
NRRD2NIFTI_BIN = "/software/DTIPrep/dev/DWI_NiftiNrrdConversion"
DWICONVERT_BIN = "/software/DTIPrep/3f0a2b4/DWIConvert"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DTIPrep on STUT data")
    parser.add_argument("sID", type=str, help="Subject ID");
    parser.add_argument("--rerun", dest="bRerun", action="store_true", \
                        help="Force rerun")
    parser.add_argument("--convertOnly", dest="bConvertOnly", \
                        action="store_true", \
                        help="Do nifti-to-nrrd conversion only")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    sID = args.sID
    bRerun = args.bRerun
    bConvertOnly = args.bConvertOnly
    
    # ------ Locate the NIFTI data ------ #
    if sID.startswith("SEQ"): # Deryk's study:
        rawdifdir = os.path.join(DATA_DIR, sID, "dwi");
        if not os.path.isdir(rawdifdir):
            raise Exception, "Cannot find raw NIFTI diffusion data directory: %s"\
                             %(rawdifdir)
        bva = glob.glob(os.path.join(rawdifdir, "*_*.mghdti.bvals"))
        if len(bva) == 0:
            raise Exception, "Found no DTI volumes in raw data directory: %s"\
                             %(rawdifdir)
        bva.sort()
        if len(bva) > 1:
            print("WARNING: found multiple DTI series in raw data directory: %s"\
                  %(rawdifdir) + " - Using the last one")
        bva = bva[-1]

        rootname = os.path.split(bva)[1].replace('.mghdti.bvals', '')
        print("INFO: rootname = %s"%(rootname))
        
        raw4d = os.path.join(rawdifdir, "%s.nii"%(rootname))
        bve = os.path.join(rawdifdir, "%s.mghdti.bvecs"%(rootname))

        if not os.path.isfile(raw4d):
            raise Exception, "Cannot find raw 4D DTI series: %s"%(raw4d)
        if not os.path.isfile(bve):
            raise Exception, "Cannot find bvecs file: %s"%(bve)

    elif len(sID) == 3 and sID.startswith("S"): # MIT STUT study
        rawdifdir = os.path.join(DATA_DIR, sID, "diffusion");
        if not os.path.isdir(rawdifdir):
            raise Exception, "Cannot find raw NIFTI diffusion data directory: %s"\
                             %(rawdifdir)

        bva = glob.glob(os.path.join(rawdifdir, "%s_*.mghdti.bvals"%(sID)))
        if len(bva) == 0:
            raise Exception, "Found no DTI volumes in raw data directory: %s"\
                             %(rawdifdir)
        bva.sort()
        if len(bva) > 1:
            print("WARNING: found multiple DTI series in raw data directory: %s"\
                  %(rawdifdir) + " - Using the last one")
        bva = bva[-1]
        rootname = os.path.split(bva)[1].replace('.mghdti.bvals', '')
        print("INFO: rootname = %s"%(rootname))

        raw4d = os.path.join(rawdifdir, "%s.nii.gz"%(rootname))
        bve = os.path.join(rawdifdir, "%s.mghdti.bvecs"%(rootname))

        if not os.path.isfile(raw4d):
            raise Exception, "Cannot find raw 4D DTI series: %s"%(raw4d)
        if not os.path.isfile(bve):
            raise Exception, "Cannot find bvecs file: %s"%(bve)
    else:
        DATA_DIR = DATA_DIR_RHY
        OUTPUT_DIR = OUTPUT_DIR_RHY
        rawdifdir = os.path.join(DATA_DIR, sID, "diffusion")

        if not os.path.isdir(rawdifdir):
            raise Exception, "Cannot find raw NIFTI diffusion data directory: %s"\
                             %(rawdifdir)

        bva = glob.glob(os.path.join(rawdifdir, "%s_diffusion_*.bval" % sID))
        if len(bva) == 0:
            raise Exception, "Found no DTI volumes in raw data directory: %s"\
                             %(rawdifdir)
        
        bva = bva[-1]
        rootname = os.path.split(bva)[1].replace('.bval', '')
        print("INFO: rootname = %s"%(rootname))

        raw4d = os.path.join(rawdifdir, "%s.nii.gz"%(rootname))
        check_file(raw4d)
        
        bve = os.path.join(rawdifdir, "%s.bvec"%(rootname))
        check_file(bve)

    print("INFO: raw4d = %s"%(raw4d))
    print("INFO: bva = %s"%(bva))
    print("INFO: bve = %s"%(bve))

    # --- Convert to nrrd --- #
    qcDir = os.path.join(OUTPUT_DIR, sID)
    if not os.path.isdir(qcDir):
        os.system("mkdir %s"%qcDir)
        print("Created directory: %s"%(qcDir))

    nrrdfn = os.path.join(qcDir, "%s.nrrd"%(sID))

    bDoCvt = bRerun or (not os.path.isfile(nrrdfn))

    if (not bRerun) and os.path.isfile(nrrdfn):
        print("INFO: nrrd file already exists: %s"%(nrrdfn))
    if bDoCvt:
        cvt_cmd = "%s --inputVolume %s "%(DWICONVERT_BIN, raw4d) + \
                  "--inputBValues %s "%(bva) + \
                  "--inputBVectors %s "%(bve) + \
                  "--conversionMode FSLToNrrd " + \
                  "--outputVolume %s"%(nrrdfn)
        saydo(cvt_cmd)
        if not os.path.isfile(nrrdfn):
            raise Exception, "Failed to convert nifti to nrrd: %s"%(nrrdfn)
        else:
            print("INFO: nifti-to-nrrd conversion done: output = %s\n"%(nrrdfn))

    if bConvertOnly:
        sys.exit(0)

    # ------- Convert to Nrrd ------ #
    # Find the diffusion_dcm directory
    '''
    rawdcmdir = os.path.join(DATA_DIR, sID, "diffusion_dcm")
    if not os.path.isdir(rawdcmdir):
        raise Exception, "Cannot find raw DTI DICOM directory: %s"%rawdcmdir
    print("INFO: rawdcmdir = %s"%rawdcmdir)
    
    qcDir = os.path.join(OUTPUT_DIR, sID)
    if not os.path.isdir(qcDir):
        os.system("mkdir %s"%qcDir)
        print("Created directory: %s"%(qcDir))

    nrrdfn = os.path.join(qcDir, "%s.nrrd"%(sID))    
    if os.path.isfile(nrrdfn) and (not bRerun):
        print("INFO: nrrd-format raw 4D data already exists: %s"%nrrdfn)
    else:
        cvt_cmd = "%s/DicomToNrrdConverter --inputDicomDirectory %s "\
                  %(DTIPREP_PATH, rawdcmdir) + \
                  "--outputDirectory %s "%(qcDir) + \
                  "--outputVolume %s.nrrd "%(sID)
        saydo(cvt_cmd)

        if not os.path.isfile(nrrdfn):
            raise Exception, "Converting dicom to nrrd failed: cannot find file: %s"%(nrrdfn)
    '''
            
    # ------ Run DTIPrep ------ #
    qced_nrrd = os.path.join(qcDir, "%s_QCed.nrrd"%(sID))
    qced_bva = os.path.join(qcDir, "%s_QCed.bvals"%(sID))
    qced_bve = os.path.join(qcDir, "%s_QCed.bvecs"%(sID))

    if os.path.isfile(qced_nrrd) and os.path.isfile(qced_bva) \
       and os.path.isfile(qced_bve) and (not bRerun):
        print("INFO: It appears that DTIprep has already been run.")
    else:
        cwd = os.getcwd()
        os.chdir(qcDir)
        dtiprep_cmd = "%s/DTIPrep -w %s.nrrd -p default -d -c "\
                      %(DTIPREP_PATH, sID) + \
                      "-n dtiprep-working.txt"
        saydo(dtiprep_cmd)
        os.chdir(cwd)

        if not os.path.isfile(qced_nrrd):
            raise Exception,  "It appears that DTIPrep has failed: cannot find series: %s"%(qced_nrrd)

    sys.exit(0)

    # ------ Convert QCed back to nifti ------ #
    qced_ngz = os.path.join(qcDir, "%s_QCed.nii.gz"%(sID))
    
    if os.path.isfile(qced_ngz) and (not bRerun):
        print("INFO: DTIPrep output has already been converted to nifti: %s"\
              %(qced_ngz))
    else:
        icvt_cmd = "%s -i %s -o %s"%(NRRD2NIFTI_BIN, qced_nrrd, qced_ngz)
        saydo(icvt_cmd)

    # ------ Run dt_recon (FreeSurfer) ------ #
    # Check the existence of the FreeSurfer data
    sFSDir = os.path.join(FSDATA_DIR, sID)
    if not os.path.isdir(sFSDir):
        raise Exception, "Cannot find FreeSurfer data directory of the subject: %s"%sFSDir
    else:
        print("INFO: sFSDir = %s"%(sFSDir))
    
    dt_recon_dir = os.path.join(OUTPUT_DIR, sID, "dt_recon")
    if not os.path.isdir(dt_recon_dir):
        os.system("mkdir %s"%dt_recon_dir)
        print("Created directory: %s"%(dt_recon_dir))

    dt_recon_cmd = "dt_recon --i %s --b %s %s --s %s --o %s"\
                   %(qced_ngz, qced_bva, qced_bve, sID, dt_recon_dir)
    saydo(dt_recon_cmd)
