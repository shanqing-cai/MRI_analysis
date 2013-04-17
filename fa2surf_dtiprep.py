#!/usr/bin/python

import os
import sys
import argparse
from scai_utils import *

FSDATA_dir = '/users/cais/STUT/FSDATA'
TRACULA_dir = "/users/cais/STUT/analysis/dti2/tracula"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Project DTIPrep'd FA data from volume to FreeSurfer surfaces")
    ap.add_argument("sID", help="Subject ID")
    ap.add_argument("projDist", help="Projection depth (Must be a negative number, e.g., -2)")
    ap.add_argument("fwhmVol", help="Volume smoothing FWHM (recommended: 0)")
    ap.add_argument("fwhmSurf", help="Surface smoothing FWHM (e.g., 5, 10, 15)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
    
    args = ap.parse_args()
    sID = args.sID
    projDist = args.projDist
    fwhmVol = args.fwhmVol
    fwhmSurf = args.fwhmSurf

    check_dir(FSDATA_dir)
    check_dir(TRACULA_dir)

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

    # Locate the fa file
    fa_file = os.path.join(subjDir, 'dmri', 'dtifit_FA.nii.gz')
    if not os.path.isfile(fa_file):
        raise IOError, "FA file not found: %s"%fa_file

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
                                "%s.FA_%dmm_fwhm%d.mgh" \
                                % (hemi, -int(projDist), int(fwhmSurf)))

        projCmd = "mri_vol2surf --mov %s --reg %s --o %s --hemi %s --trgsubject %s --projdist %s --fwhm %s --surf-fwhm %s --noreshape --interp trilin" \
            % (fa_file, reg_dat, out_file, hemi, \
               sID, projDist, fwhmVol, fwhmSurf)

        os.system("rm -f %s" % out_file)

        saydo(projCmd)
        check_file(out_file)

        viewCmds[hemi] = \
            'tksurfer %s %s inflated -gray -overlay %s -fminmax 0.01 1' \
            % (sID, hemi, out_file)

        fsaverage_fn = os.path.join(FSSubjDir, 'surf', \
                                    '%s.FA_%dmm_fwhm%d.fsaverage.mgh'\
                                    %(hemi, -int(projDist), int(fwhmSurf)))

        preprocCmd = "mris_preproc --s %s --hemi %s --meas FA_%dmm_fwhm%d.mgh --target fsaverage --out %s" \
                     %(sID, hemi, -int(projDist), int(fwhmSurf), fsaverage_fn)
        
        os.system("rm -f %s" % fsaverage_fn)

        saydo(preprocCmd)
        check_file(fsaverage_fn)
        
    print('----------------------------------------------------------')
    for hemi in hemis:
        print('To view the result of %s: \n\t%s\n'%(hemi, viewCmds[hemi]))

