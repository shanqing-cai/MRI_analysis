#!/usr/bin/python

import os
import sys
import argparse
from scai_utils import *

FSDATA_dir = '/users/cais/STUT/FSDATA'
TRACULA_dir = "/users/cais/STUT/analysis/dti2/tracula"

VALID_MEAS = ["FA", "L1", "RD"];

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Project DTIPrep'd FA data from volume to FreeSurfer surfaces")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("projDist", help="Projection depth (Must be a negative number, e.g., -2)")
    ap.add_argument("fwhmVol", help="Volume smoothing FWHM (recommended: 0)")
    ap.add_argument("fwhmSurf", help="Surface smoothing FWHM (e.g., 5, 10, 15)")
    ap.add_argument("--meas", dest="meas", default="FA", 
                    help="Diffusion tensor measure {L1, RD} (RD is defined as the arithmetic mean of L2 and L3)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
    
    args = ap.parse_args()
    sID = args.sID
    projDist = args.projDist
    meas =args.meas
    fwhmVol = args.fwhmVol
    fwhmSurf = args.fwhmSurf

    check_dir(FSDATA_dir)
    check_dir(TRACULA_dir)

    # Input arguments sanity check
    if VALID_MEAS.count(meas) == 0:
        raise Exception, "Unrecognized diffusion tensor measure: %s" % meas

    # Locate the FS subject
    FSSubjDir = os.path.join(FSDATA_dir, sID)
    if not os.path.isdir(FSSubjDir):
        raise IOError, "Subject %s not found in FreeSurfer directory %s"\
            %(sID, FSDATA_dir)

    # Locate the dt_recon subject
    subjDir = os.path.join(TRACULA_dir, sID)
    if not os.path.isdir(subjDir):
        raise IOError, "Subject %s not found in dt_recon directory %s"\
            %(sID, dt_recon_dir)

    # Locate the diffusion-tensor measure file
    if meas == "FA" or meas == "L1":
        dtmeas_file = os.path.join(subjDir, 'dmri', \
                                   'dtifit_%s.nii.gz' % meas)
        if not os.path.isfile(dtmeas_file):
            raise IOError, "Diffusion-tensor measure (%s) file not found: %s" \
                           % (meas, dtmeas_file)
    elif meas == "RD":
        dtmeas_file = os.path.join(subjDir, "dmri", \
                                   "dtifit_RD.nii.gz")
        L2_file = os.path.join(subjDir, "dmri", "dtifit_L2.nii.gz")
        L3_file = os.path.join(subjDir, "dmri", "dtifit_L3.nii.gz")
        check_file(L2_file)
        check_file(L3_file)

        t_sum = os.path.join(subjDir, "dmri", "dtifit_RDt2.nii.gz")
        saydo("rm -f %s" % t_sum)
        add_cmd = "fslmaths %s -add %s %s" % \
                  (L2_file, L3_file, t_sum)
        saydo(add_cmd)
        check_file(t_sum)

        saydo("rm --f %s" % dtmeas_file)
        div_cmd = "fslmaths %s -div 2.0 %s" % (t_sum, dtmeas_file)
        saydo(div_cmd)
        check_file(dtmeas_file)


    # Locate the diff2anat reg data
    reg_dat = os.path.join(subjDir, "dmri", "xfms", "d2a.dat")
    if not os.path.isfile(reg_dat):
        raise IOError, \
              "register.dat not found: %s" % reg_dat

    if float(projDist) >= 0:
        raise IOError, "projDist should be a negative number"

    if float(fwhmVol) < 0:
        raise IOError, "fwhmVol should be zero or positive"

    if float(fwhmSurf) < 0:
        raise IOError, "fwhmSurf should be zero or positive"

    hemis = ['lh', 'rh']
    viewCmds = {'lh': '', 'rh': ''}

    for hemi in hemis:
        out_file = os.path.join(FSSubjDir, "surf", \
                                "%s.%s_%dmm_fwhm%d.mgh" \
                                % (hemi, meas, -int(projDist), int(fwhmSurf)))

        projCmd = "mri_vol2surf --mov %s --reg %s --o %s --hemi %s --trgsubject %s --projdist %s --fwhm %s --surf-fwhm %s --noreshape --interp trilin" \
            % (dtmeas_file, reg_dat, out_file, hemi, \
               sID, projDist, fwhmVol, fwhmSurf)

        os.system("rm -f %s" % out_file)

        saydo(projCmd)
        check_file(out_file)

        viewCmds[hemi] = \
            'tksurfer %s %s inflated -gray -overlay %s -fminmax 0.01 1' \
            % (sID, hemi, out_file)

        fsaverage_fn = os.path.join(FSSubjDir, 'surf', \
                                    '%s.%s_%dmm_fwhm%d.fsaverage.mgh'\
                                    %(hemi, meas, \
                                      -int(projDist), int(fwhmSurf)))

        preprocCmd = "mris_preproc --s %s --hemi %s --meas FA_%dmm_fwhm%d.mgh --target fsaverage --out %s" \
                     %(sID, hemi, -int(projDist), int(fwhmSurf), fsaverage_fn)
        
        os.system("rm -f %s" % fsaverage_fn)

        saydo(preprocCmd)
        check_file(fsaverage_fn)
        
    print('----------------------------------------------------------')
    for hemi in hemis:
        print('To view the result of %s: \n\t%s\n'%(hemi, viewCmds[hemi]))

