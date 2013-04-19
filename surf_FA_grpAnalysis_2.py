#!/usr/bin/python

import os
import sys
import glob
import argparse
from get_qdec_info import *
from scai_utils import *
from exclusion_paradigms import exclParad

# Usage: surf_FA_grpAnalysis.py fwhm subjsList
FSDATA_dir = '/users/cais/STUT/FSDATA/'
qdec_dir = '/users/cais/STUT/FSDATA/'
surf_FA_dir = '/users/cais/STUT/analysis/surf_FA_grp'
MATLAB_BIN = '/software/matlab2007b/bin/matlab'

con_10 = '/users/cais/STUT/analysis/design/con_10'
con_01 = '/users/cais/STUT/analysis/design/con_01'

SURF_CSD_WC = "/software/Freesurfer/5.0.0/average/mult-comp-cor/fsaverage/rh/cortex/%s/%s/%s/mc-z.csd"

VALID_MEAS = ["FA", "L1", "RD"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Group-level analysis")
    ap.add_argument("subjsList", \
        help="Subject list (e.g., /users/cais/STUT/FSDATA/qdec/subjsList.txt")
    ap.add_argument("depth", \
        help="Projection depth (e.g., 2)")
    ap.add_argument("fwhm", \
        help="Surface FWHM (e.g., 15)")
    ap.add_argument("--meas", dest="meas", default="FA", \
        help="Type of diffusion-tensor measure {FA, L1, RD}")
    ap.add_argument("--exclParad", dest="exclName", default="", \
        help="Exclusion paradigm (e.g., noS24)")
    ap.add_argument("--clusterP", type=str, dest="clustPCfg", default="", \
        help="Estimate cluster p-values. Requires input: surface FWHM, vertex threshold, sign, cluster-wise p-value threshold (e.g., 5,0.01,abs,0.05)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    subjsList = args.subjsList
    depth = args.depth
    fwhm = args.fwhm
    meas = args.meas
    exclName = args.exclName
    clustPCfg = args.clustPCfg

    # Input sanity check
    if VALID_MEAS.count(meas) == 0:
        raise Exception, "Unrecognized diffusion tensor measure: %s" % meas

    # --- If cluster p is requested --- #
    # --- Locate the proper csd file --- #
    if len(clustPCfg) > 0:
        if clustPCfg.count(',') != 3:
            raise Exception, "Erroneous syntax for clustPCfg. Correc tusage: surfFWHM,vThr,sign (e.g., 5,0.01,abs)"

        surfFWHM = int(clustPCfg.split(',')[0])
        if surfFWHM < 0:
            raise Exception, "clustPCfg: surfFWHM must be a non-negative number"

        if surfFWHM != int(fwhm):
            raise Exception, "The value of clustPCfg:surfFWHM (%d) does not match the fwhm in main input arguments (%d)" \
                             % (surfFWHM, int(fwhm))
        
        vThr = float(clustPCfg.split(',')[1])
        if vThr <= 0 or vThr >= 1:
            raise Exception, "clustPCfg: surfFWHM value error (%f)" % vThr
        
        clustSign = clustPCfg.split(',')[2]
        
        fwhm_surf_str = "fwhm%.2d" % (surfFWHM)
        
        cwpThresh = float(clustPCfg.split(',')[3])
        if cwpThresh <= 0.0 or cwpThresh >= 1.0:
            raise Exception, "clustPCfg: cwpThresh value error (%f)" % cwpThresh
        
        if vThr == 0.05:
            th_str = "th13"
        elif vThr == 0.01:
            th_str = "th20"
        elif vThr == 0.005:
            th_str = "th23"
        elif vThr == 0.001:
            th_str = "th30"
        elif vThr == 0.0005:
            th_str = "th33"
        elif vThr == 0.00001:
            th_str = "th40"
        else:
            raise Exception, "Unrecognized value of MCZ_SURF_VERTEX_THRESH: %f" % vThr

        csdfn = SURF_CSD_WC % (fwhm_surf_str, clustSign, th_str)
        check_file(csdfn)


    f_subjsList = open(subjsList, 'r')
    t_subjs = f_subjsList.read()
    f_subjsList.close()

    subjs = t_subjs.split('\n')
    subjs.remove('')
    
    #--- Apply exclusion paradigm ---#
    if exclName != '':
        exclS = exclParad[exclName]
        for e_sID in exclS:
            if subjs.count(e_sID) == 1:
                subjs.remove(e_sID)
                print('Removing subject %s from this analysis'%(e_sID))
    #--- ~Apply exclusion paradigm ---#    

    nSubjs = len(subjs)

    diags = [''] * nSubjs
    diagVec = [0] * nSubjs
    X_line = 'X = ['

    if meas == "L1":
        surf_FA_dir = surf_FA_dir.replace("_FA_", "_L1_")
    elif meas == "RD":
        surf_FA_dir = surf_FA_dir.replace("_FA_", "_RD_")

    if not os.path.isdir(surf_FA_dir):
        raise IOError, 'Directory not found: %s'%surf_FA_dir
    X_mat = os.path.join(surf_FA_dir, 'X.mat');

    for i1 in range(nSubjs):
        s = subjs[i1]

        if len(s) == 3:
            diags[i1] = get_qdec_info(s, 'diagnosis')
        elif len(s) == 8:
            if s[5] == 'C':
                diags[i1] = "PFS"
            elif s[5] == 'P':
                diags[i1] = "PWS"
            else:
                raise Exception, "Unrecognized subject ID format: %s"%s
        else:
            raise Exception, "Unrecognized subject ID format: %s"%s

        diagVec[i1] = int(((diags[i1] == 'PWS') - 0.5) * 2)
        # So +1 for PWS and -1 for PFS

        X_line += '%d; '%diagVec[i1]

    if diags.count(''):
        raise Exception, "Failed to determine the group identity of %d subjects"%(diags.count(''))

    # Generate the 
    X_line = X_line[:-2]
    X_line += '];\nX = [ones(' + str(nSubjs) + ', 1), X];'

    X_mat = os.path.join(surf_FA_dir, 'X_%s.mat'%exclName)
    X_line += '\nsave(\'' + X_mat + '\', \'X\', \'-v4\')\n'

    mScriptGenX = os.path.join(surf_FA_dir, 'genX.m')
    genX_f = open(mScriptGenX, 'w')
    genX_f.write(X_line)
    genX_f.close()

    matlabCmd = MATLAB_BIN + ' -nosplash -nodesktop -r \'run ' \
                           + mScriptGenX + '; exit\''
    saydo(matlabCmd)
    check_file(X_mat)

    # Generate the merged surf_FA files
    hemis = ['lh', 'rh']
    viewCmds = {'lh': '', 'rh': ''}

    cwsig_fns = {"lh": "", "rh": ""}
    vwsig_fns = {"lh": "", "rh": ""}
    ocn_fns = {"lh": "", "rh": ""}
    mcz_surf_sum_fns = {"lh": "", "rh": ""}

    for hemi in hemis:
        concatCmd = 'mri_concat --i '
        for s in subjs:
            surf_FA_fn = os.path.join(FSDATA_dir, s, 'surf', \
                                      '%s.%s_%smm_fwhm%s.fsaverage.mgh'\
                                      %(hemi, meas, depth, fwhm))

            if not os.path.isfile(surf_FA_fn):
                raise IOError, 'Data file not found: %s'%surf_FA_fn
            concatCmd += surf_FA_fn + ' '

        mergedSurfFA = os.path.join(surf_FA_dir, \
                                    '%s.merged_%s_%dS_%smm_fwhm%s.fsaverage.mgh'
                                    %(hemi, meas, nSubjs, depth, fwhm))

        concatCmd += '--o ' + mergedSurfFA
        saydo(concatCmd)

        if not os.path.isfile(mergedSurfFA):
            raise Exception, "It appears that mri_concat has failed to generate the merged FA file for %s: %s"%(hemi, mergedSurfFA)

        if exclName == '':
            con_10_dir = os.path.join(surf_FA_dir, 'con_10_%smm_fwhm%s.%s'\
                                                   %(depth, fwhm, hemi))
        else:
            con_10_dir = os.path.join(surf_FA_dir, 'con_10_%smm_fwhm%s_%s.%s'\
                                                   %(depth, fwhm, exclName, hemi))        
    
        fitCmd = 'mri_glmfit --y ' + mergedSurfFA + ' --X ' + X_mat + ' --C ' \
                 + con_10 + ' --glmdir ' + con_10_dir
        saydo(fitCmd)
        check_dir(con_10_dir)
        check_file(os.path.join(con_10_dir, "con_10", "sig.mgh"))

        if exclName == '':
            con_01_dir = os.path.join(surf_FA_dir, 'con_01_%smm_fwhm%s.%s'\
                                                   %(depth, fwhm, hemi))
        else:
            con_01_dir = os.path.join(surf_FA_dir, 'con_01_%smm_fwhm%s_%s.%s'\
                                                   %(depth, fwhm, exclName, hemi))        
            
        fitCmd = 'mri_glmfit --y ' + mergedSurfFA + ' --X ' + X_mat + ' --C ' \
                 + con_01 + ' --glmdir ' + con_01_dir

        saydo(fitCmd)
        check_dir(con_01_dir)
        sigFile = os.path.join(con_01_dir, "con_01", "sig.mgh")
        check_file(sigFile)    

        viewCmds[hemi] = 'tksurfer fsaverage %s inflated -gray -overlay %s -fminmax 3 6'%(hemi, sigFile)

        # --- Optional: surface cluster p-value calculation --- #
        # --- (Using FreeSurfer pre-generated csd files --- #
        if len(clustPCfg) > 0:
            cwsig_fns[hemi] = os.path.join(con_01_dir, \
                                    "mc-z.abs.%s.%s.sig.cluster.mgh" % (fwhm_surf_str, th_str))
            vwsig_fns[hemi] = os.path.join(con_01_dir, \
                                               "mc-z.abs.%s.%s.sig.vertex.mgh" % (fwhm_surf_str, th_str))
            mcz_surf_sum_fns[hemi] = os.path.join(con_01_dir, \
                                              "mc-z.abs.%s.%s.cluster.sum"\
                                              % (fwhm_surf_str, th_str))
            ocn_fns[hemi] = os.path.join(con_01_dir, \
                                         "mc-z.abs.%s.%s.sig.ocn.mgh" % (fwhm_surf_str, th_str))

            surfclust_cmd = "mri_surfcluster --in %s --csd %s --cwsig %s --vwsig %s --sum %s --ocn %s --annot aparc --cwpvalthresh %f --surf white"\
                        %(sigFile, csdfn, cwsig_fns[hemi], vwsig_fns[hemi], \
                          mcz_surf_sum_fns[hemi], ocn_fns[hemi], cwpThresh)
            
            saydo(surfclust_cmd)
            check_file(cwsig_fns[hemi])
            check_file(vwsig_fns[hemi])
            check_file(mcz_surf_sum_fns[hemi])
            check_file(ocn_fns[hemi])



    print('\n================== To view the commands ==================')
    for hemi in hemis:
        if hemi == 'lh':
            sideWord = 'Left'
        else:
            sideWord = 'Right'

        print('%s: hemisphere:\n\t %s\n'%(sideWord, viewCmds[hemi]))

    print('==========================================================')

    if len(clustPCfg) > 0:
        print('\n================== Surface cluster summary files  ==================')
        for hemi in hemis:
            if hemi == 'lh':
                sideWord = 'Left'
            else:
                sideWord = 'Right'

            print("%s hemisphere:\n\tless %s\n" \
                  % (sideWord, mcz_surf_sum_fns[hemi]))


