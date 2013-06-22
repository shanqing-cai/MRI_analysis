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
MNI152_TEMPLATE_FN = "/usr/share/fsl/data/standard/MNI152_T1_2mm.nii.gz"
L2_DIR = "/users/cais/STUT/analysis/tract_subcort_conn"
DESIGN_DIR = "/users/cais/STUT/analysis/design"

MATLAB_BIN = '/software/matlab2009a/bin/matlab'
con_01 = '/users/cais/STUT/analysis/design/con_01'

regions = ['Prefrontal', 'Premotor', 'Precentral', 'Postcentral', \
           'PPC', 'Occipital', 'Temporal']

mainMaskRatio = 0.75

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Group-level analysis of subcortical-cortical DTT connectivity")
    ap.add_argument("scSeed", help="Subcortical seed")
    ap.add_argument("--ccstop", dest="bCCStop", \
                    action="store_true", help="Use corpus callosum stop mask")
    ap.add_argument("--roi", dest="roi", type=str, default="", \
                    help="Focus on specific cortical ROI (default: False, i.e., use all %d regions" % len(regions))
    ap.add_argument("--rerun", dest="bRerun", \
                    action="store_true", help="Rerun time-consuming steps")
    ap.add_argument("--mask", dest="mask", type=str, default="", \
                    help="Mask for mri_glmfit")
    ap.add_argument("--fwhm", dest="fwhm", type=float, default=0, \
                    help="FWHM of smoothing for mri_glmfit")
    ap.add_argument("--mainMaskRatio", dest="mainMaskRatio", type=float, \
                    default=mainMaskRatio, \
                    help="Ratio for generating the main PFS mask (default=%f)" \
                         % mainMaskRatio)
    
    # === Process input arguments === #
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    
    scSeed = args.scSeed
    bCCStop = args.bCCStop
    bRerun = args.bRerun
    roi = args.roi
    mask = args.mask
    fwhm = args.fwhm
    mainMaskRatio = args.mainMaskRatio

    assert(mainMaskRatio >= 0.0 and mainMaskRatio < 1.0)

    if roi == "":
        a_regions = regions
    else:
        a_regions = [roi]

    assert(scSeed.count(".") == 1)
    hemi = scSeed.split(".")[0]
    cSeed = scSeed.split(".")[1]

    if mask != "":
        check_file(mask)

    if fwhm != 0:
        assert(fwhm > 0)

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
    for i0 in range(len(a_regions)):
        regMean_ts.append([])

    for (i0, t_sID) in enumerate(sIDs):
        sDir = os.path.join(tractSegDir, t_sID, scDirName)
        check_dir(sDir)

        for (i1, t_region) in enumerate(a_regions):
            if roi == "":
                regMean_d = os.path.join(sDir, \
                                         "tract_regionMean_%d.nii.gz" % (i1))
                regMean_t = os.path.join(sDir, \
                                     "tract_regionMean_%d.MNI152.nii.gz" % (i1))
            else:
                regMean_d = os.path.join(sDir, "aparc12", 
                                         "%s.nii.gz" % roi)
                regMean_t = os.path.join(sDir, "aparc12",
                                     "%s.MNI152.nii.gz" % roi)
            check_file(regMean_d)
            
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
    regMean_4ds_PWS = []
    regMean_4ds_PFS = []
    mean_PWS = []
    mean_PFS = []
    mainMask_PFS = []
    L2SeedDir = os.path.join(L2_DIR, scDirName)
    check_dir(L2SeedDir, bCreate=True)

    for (i0, t_region) in enumerate(a_regions):
        if roi == "":
            L2SeedTargDir = os.path.join(L2SeedDir, "%d_%s" % (i0, t_region))
        else:
            L2SeedTargDir = os.path.join(L2SeedDir, t_region)
        check_dir(L2SeedTargDir, bCreate=True)
            
        if roi == "":
            regMean_4ds.append(os.path.join(L2SeedTargDir, \
                                            "merged_%d.nii.gz" % i0))
            regMean_4ds_PWS.append(os.path.join(L2SeedTargDir, \
                                            "merged_%d_PWS.nii.gz" % i0))
            regMean_4ds_PFS.append(os.path.join(L2SeedTargDir, \
                                            "merged_%d_PFS.nii.gz" % i0))
        else:
            regMean_4ds.append(os.path.join(L2SeedTargDir, \
                                            "merged_%s.nii.gz" % roi))
            regMean_4ds_PWS.append(os.path.join(L2SeedTargDir, \
                                            "merged_%s_PWS.nii.gz" % roi))
            regMean_4ds_PFS.append(os.path.join(L2SeedTargDir, \
                                            "merged_%s_PFS.nii.gz" % roi))

        mergeCmd = "fslmerge -t %s " % regMean_4ds[-1]
        mergeCmd_PWS = "fslmerge -t %s " % regMean_4ds_PWS[-1]
        mergeCmd_PFS = "fslmerge -t %s " % regMean_4ds_PFS[-1]

        for (i1, t_fn) in enumerate(regMean_ts[i0]):
            mergeCmd += "%s " % t_fn
            if bPWS[i1] == 1:
                mergeCmd_PWS += "%s " % t_fn
            else:
                mergeCmd_PFS += "%s " % t_fn
        
        mean_PWS.append(os.path.join(L2SeedTargDir, \
                                     "merged_%s_mean_PWS.nii.gz" % roi))
        mean_PFS.append(os.path.join(L2SeedTargDir, \
                                     "merged_%s_mean_PFS.nii.gz" % roi))

        mainMask_PFS.append(os.path.join(L2SeedTargDir, \
                                     "mainMask_%s_PFS_%.2f.nii.gz" \
                                     % (roi, mainMaskRatio)))
        if bRerun or (not os.path.isfile(regMean_4ds[-1]) \
                      or not os.path.isfile(regMean_4ds_PWS[-1]) \
                      or not os.path.isfile(regMean_4ds_PFS[-1]) \
                      or not os.path.isfile(mean_PWS[-1]) \
                      or not os.path.isfile(mean_PFS[-1]) \
                      or not os.path.isfile(mainMask_PFS[-1])):
            saydo(mergeCmd)
            check_file(regMean_4ds[-1])

            saydo(mergeCmd_PWS)
            check_file(regMean_4ds_PWS[-1])

            saydo(mergeCmd_PFS)
            check_file(regMean_4ds_PFS[-1])
            
            # Calculate the mean:
            meanCmd = "fslmaths %s -Tmean %s" % \
                      (regMean_4ds_PWS[-1], mean_PWS[-1])
            saydo(meanCmd)
            check_file(mean_PWS[-1])

            meanCmd = "fslmaths %s -Tmean %s" % \
                      (regMean_4ds_PFS[-1], mean_PFS[-1])
            saydo(meanCmd)
            check_file(mean_PFS[-1])

            # Get the 99 percentile of mean_PFS
            (so, se) = cmd_stdout("fslstats %s -P 99" % (mean_PFS[-1]))
            assert(len(se) == 0)
            p99 = float(so.split(' ')[0])
            print("INFO: PFS 99 percentile = %f" % p99)

            # Generate the main mask
            binCmd = "mri_binarize --i %s --min %f --binval 1.0 --o %s" \
                     % (mean_PFS[-1], p99 * mainMaskRatio, mainMask_PFS[-1])
            saydo(binCmd)
            check_file(mainMask_PFS[-1])
            
        # == Prepare for stats on the main-masked mean == #
        (so, se) = cmd_stdout("fslstats -t %s -k %s -M" % \
                              (regMean_4ds[-1], mainMask_PFS[-1]))
        #assert(len(se) == 0)
        mainMaskVals = np.zeros(len(bPWS))
        so = remove_empty_strings(so.replace("\n", " ").split(" "))
        for (i1, s) in enumerate(so):
            mainMaskVals[i1] = float(s)

        v_PWS = mainMaskVals[np.nonzero(bPWS == 1)]
        v_PFS = mainMaskVals[np.nonzero(bPWS == 0)]

        # == Run stats == #
        import scipy.stats as stats
        (tt_t, tt_p) = stats.ttest_ind(v_PWS, v_PFS)
        (rs_z, rs_p) = stats.ranksums(v_PWS, v_PFS)

        stfn = os.path.join(L2SeedTargDir, "stats_%.2f.txt" % mainMaskRatio)
        stf = open(stfn, "wt")
        stf.write("Main-masked values: \n")
        stf.write("PWS: mean=%f; ste=%f\n" \
                  % (np.mean(v_PWS), np.std(v_PWS) / np.sqrt(len(v_PWS))))
        for (i1, v) in enumerate(v_PWS):
            stf.write("%f " % v)
        stf.write("\n")

        stf.write("PFS: mean=%f; ste=%f\n" \
                  % (np.mean(v_PFS), np.std(v_PFS) / np.sqrt(len(v_PFS))))
        for (i1, v) in enumerate(v_PFS):
            stf.write("%f " % v)
        stf.write("\n\n")
        
        stf.write("t-test:\n\tt=%f; p=%f\n" % (tt_t, tt_p))
        stf.write("ranksum:\n\tZ=%f; p=%f\n" % (rs_z, rs_p))
        
        stf.close()

        check_file(stfn)
        print("INFO: statistical results written to file:\n\t%s" % stfn)
        

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
    
    for (i0, t_region) in enumerate(a_regions):
        if roi == "":
            t_con_01_dir = os.path.join(L2SeedTargDir, "bgc")
        else:
            t_con_01_dir = os.path.join(L2SeedTargDir, "bgc")

        fitCmd = 'mri_glmfit --no-prune --y ' + regMean_4ds[i0] \
                 + ' --X ' + X_mat \
                 + ' --C ' + con_01 \
                 + ' --glmdir ' + t_con_01_dir

        if mask != "":
            fitCmd += " --mask %s" % mask

        if fwhm > 0:
            fitCmd += " --fwhm %f" % fwhm

        saydo(fitCmd)
        check_dir(t_con_01_dir)

        sig_fn = os.path.join(t_con_01_dir, "con_01", "sig.mgh")
        check_file(sig_fn)
    
        # == Convert the sig file to nii.gz == #
        sig_ngz = sig_fn.replace(".mgh", ".nii.gz")
        cvtCmd = "mri_convert %s %s" % \
                 (sig_fn, sig_ngz)
        saydo(cvtCmd)
        check_file(sig_ngz)

        
        
