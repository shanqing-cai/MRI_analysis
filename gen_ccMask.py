# gen_ccMask.py: Generate corpus callosum mask from FreeSufer autoseg
# Author: Shanqing Cai
# Date: 2013-02-17

import os
import sys

if not len(sys.argv) == 2:
    print('Usage: gen_ccMask.py sID')
    sys.exit(0)

FSBASE = '/users/cais/STUT/FSDATA'
TRACBASE = '/users/cais/STUT/analysis/dti2/tracula'

sID = sys.argv[1]

slaparc_fn = os.path.join(FSBASE, sID, 'mri', 'slaparc.nii.gz')

if not os.path.isfile(slaparc_fn):
    raise IOError, 'Cannot find slaparc file %s'%slaparc_fn

print('slaparc_fn = %s'%slaparc_fn)

ccMask_fn = os.path.join(FSBASE, sID, 'mri', 'ccMask.nii.gz')
binarize_cmd = 'mri_binarize --i %s --min 251 --max 255 --o %s'\
               %(slaparc_fn, ccMask_fn)

def saydo(cmd):
    print(cmd)
    os.system(cmd)

saydo(binarize_cmd)

anat2diff_fsl = os.path.join(TRACBASE, sID, 'dmri', 'xfms', \
                                 'anatorig2diff.bbr.mat')
if not os.path.isfile(anat2diff_fsl):
    raise IOError, 'Cannot find anat2diff FSL xfm: %s'%anat2diff_fsl

ref_img = os.path.join(TRACBASE, sID, 'dmri', 'dtifit_FA.nii.gz')
if not os.path.isfile(anat2diff_fsl):
    raise IOError, 'Cannot find FA image: %s'%ref_img

ccMask_diff_fn = os.path.join(FSBASE, sID, 'mri', 'ccMask_diff.nii.gz')
xfm_cmd = 'flirt -in %s -ref %s -applyxfm -init %s -out %s -interp trilinear'\
          %(ccMask_fn, ref_img, anat2diff_fsl, ccMask_diff_fn)

saydo(xfm_cmd)


rebin_cmd = 'mri_binarize --i %s --min 0.001 --max 1 --o %s'\
            %(ccMask_diff_fn, ccMask_diff_fn)
saydo(rebin_cmd)

print('\nccMask in structural space: %s'%ccMask_fn)
print('ccMask in diffusion space: %s\n'%ccMask_diff_fn)
