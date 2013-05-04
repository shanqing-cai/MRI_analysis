#!/usr/bin/python

import sys
import os
import glob
import argparse
import numpy as np

from scai_utils import *
from get_qdec_info import get_qdec_info

## CONFIG: Paths
#qdec_dat = '/users/cais/STUT/FSDATA/qdec/qdec.table.dat'
qdec_comp = '/users/cais/STUT/FSDATA/qdec/qdec_PWS.dat'

thresh_SSI4 = 16
SKEL_FA_THRESH = 0.2

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="TBSS statistical analysis")
    ap.add_argument("tbss_dir", type=str, \
                    help="Base TBSS directory (e.g., /users/cais/STUT/analysis/dt_tbss_dtiprep2")
    ap.add_argument("tbss_dir_new", type=str, \
                    help="Output TBSS directory")
    ap.add_argument("--incGen", dest="bIncGen", action="store_true", \
                    help="Include gender as a covariate of non-interest")
    ap.add_argument("--incAge", dest="bIncAge", action="store_true", \
                    help="Include age as covariate of non-interest")
    ap.add_argument("--maleOnly", dest="bMaleOnly", action="store_true", \
                    help="Include only male subjects (exclude female subjects)")
    ap.add_argument("--noMild", dest="bNoMild", action="store_true", \
                    help="Exclude mild PWS")
    ap.add_argument("--corrSSI", dest="bCorrSSI", action="store_true", \
                    help="Perform linear correlation with SSI-4 severity score")

## Get the list of subjects
    if len(sys.argv) == 1:
#        sys.exit('Wrong syntax')
        ap.print_help()
        sys.exit(0)

    # === Parse input arguments ===
    args = ap.parse_args()
    tbss_dir = args.tbss_dir
    tbss_dir_new = args.tbss_dir_new
    bMaleOnly = args.bMaleOnly
    bIncGen = args.bIncGen
    bIncAge = args.bIncAge
    bNoMild = args.bNoMild
    bCorrSSI = args.bCorrSSI

    assert(tbss_dir != tbss_dir_new)

    if bCorrSSI and bNoMild:
        raise Exception, \
              "Do NOT use the options -noMild and -corrSSI simultaneously."

    if bIncGen and bMaleOnly:
        raise Exception, \
              "Do NOT use the options -incGen and -maleOnly simultaneously."

    if bIncGen and bIncAge:
        raise Exception, \
              "Simultaneous use of --incGen and --incAge is not supported in this version"

#if noMild == '-noMild':
#    bNoMild = True
#else:
#    bNoMild = False

    if bIncAge:
        tbss_mat = '/users/cais/STUT/analysis/design/design_tbss_bgc_incAge.mat'
        tbss_con = '/users/cais/STUT/analysis/design/design_tbss_bgc_incAge.con'
        tbss_con_rev = '/users/cais/STUT/analysis/design/design_tbss_bgc_incAge_rev.con'
    elif bIncGen:
        tbss_mat = '/users/cais/STUT/analysis/design/design_tbss_bgc_incGen.mat'
        tbss_con = '/users/cais/STUT/analysis/design/design_tbss_bgc_incGen.con'
        tbss_con_rev = '/users/cais/STUT/analysis/design/design_tbss_bgc_incGen_rev.con'
    elif bMaleOnly:
        tbss_mat = '/users/cais/STUT/analysis/design/design_tbss_bgc_maleOnly.mat'
        tbss_con = '/users/cais/STUT/analysis/design/design_tbss_bgc_maleOnly.con'
        tbss_con_rev = '/users/cais/STUT/analysis/design/design_tbss_bgc_maleOnly_rev.con'
    else:
        tbss_mat = '/users/cais/STUT/analysis/design/design_tbss_bgc.mat'
        tbss_con = '/users/cais/STUT/analysis/design/design_tbss_bgc.con'
        tbss_con_rev = '/users/cais/STUT/analysis/design/design_tbss_bgc_rev.con'

    
    check_dir(tbss_dir)
#    if os.path.isdir(tbss_dir):
    print('tbss_dir = ' + tbss_dir)
#    else:
#        sys.exit('ERROR: tbss_dir ' + tbss_dir + ' not found.')

    orig_dir = os.path.join(tbss_dir, 'origdata')
    check_dir(orig_dir)

#if not os.path.isfile(qdec_comp):
#    sys.exit('ERROR: qdec_comp file ' + qdec_comp + ' not found')

    d_orig = glob.glob(os.path.join(orig_dir, 'S*.nii.gz'))
    S_list = []
    for fn in d_orig:
        S_list.append(fn.split('/')[-1].split('.')[0])
    S_list.sort()
    print('\nFound ' + str(len(S_list)) + ' subjects in ' + orig_dir + '\n')
    # === Get subjects info ===
    isPWS_list = [0] * len(S_list)
    isFemale_list = [0] * len(S_list)
    age_list = [-1.0] * len(S_list)
#    done_list = [0] * len(S_list)
#    SSI4_list = [-1.0] * len(S_list)

    for (i0, sID) in enumerate(S_list):
        t_diag = get_qdec_info(sID, "diagnosis")
        if t_diag == "PWS":
            isPWS_list[i0] = 1
        else:
            isPWS_list[i0] = 0

        t_gen = get_qdec_info(sID, "gender")
        if t_gen.lower() == "female":
            isFemale_list[i0] = 1
        else:
            isFemale_list[i0] = 0
            
        age_list[i0] = int(get_qdec_info(sID, "Age"))

    isPWS_list = np.array(isPWS_list)
    isFemale_list = np.array(isFemale_list)
    age_list = np.array(age_list)

    N = len(S_list)
    assert(len(isPWS_list) == N)
    assert(len(isFemale_list) == N)
    assert(len(age_list) == N)

    # == Center age vector ==
    age_list = age_list - np.mean(age_list)

    # === Figure out the SSI4 scores for the PWS ===
    qcf = open(qdec_comp, 'r')
    qct = qcf.read()
    qct = qct.split('\n')
    qcf.close()

    qct = qct[1:]

    if bNoMild or bCorrSSI:
        for i1 in range(len(S_list)):
            if isPWS_list[i1] == 0:
                continue
        #print('Figuring out SSI4 score for ' + S_list[i1])
            bFound = False
            for t_line in qct:
                t_line_1 = t_line.replace('\t', ' ').split(' ')

                while t_line_1.count('') > 0:
                    t_line_1.remove('')

                if len(t_line_1) < 4:
                    continue
                
                t_sID = t_line_1[0]
                if t_sID == S_list[i1] or \
                   t_sID.replace('_500', '') == S_list[i1]:
                    bFound = True
                    SSI4_list[i1] = float(t_line_1[3])
                    
                    break

            if not bFound:
                sys.exit('ERROR: SSI4 score not found for subject ' \
                         + S_list[i1])
            else:
                print(S_list[i1] + ': SSI4 = ' + str(SSI4_list[i1]))

    # === === #
    if bNoMild:         # Not including mild PWS subjects
        S_list_new = []
        isPWS_list_new = []
        isFemale_list_new = []
        S_toDel = []

        for i1 in range(len(SSI4_list)):
            if SSI4_list[i1] >= thresh_SSI4 or SSI4_list[i1] == -1:
                S_list_new.append(S_list[i1])
                isPWS_list_new.append(isPWS_list[i1])
                isFemale_list_new.append(isFemale_list[i1])
            else:
                S_toDel.append(S_list[i1])
    
    elif bMaleOnly:     ## Male only
        S_list_new = []
        isPWS_list_new = []
        isFemale_list_new = []
        S_toDel = []

        for i1 in range(len(SSI4_list)):
            if isFemale_list[i1] == 0:
                S_list_new.append(S_list[i1])
                isPWS_list_new.append(isPWS_list[i1])
                isFemale_list_new.append(isFemale_list[i1])
            else:
                S_toDel.append(S_list[i1])

    elif bCorrSSI:      ## Include SSI4 scores as a covariate
        S_list_new = [];
        isPWS_list_new = [];
        isFemale_list_new = [];
        S_toDel = [];
        SSI_reg = [];

        for i1 in range(len(SSI4_list)):
            if SSI4_list[i1] >= 0:
                S_list_new.append(S_list[i1])
                isPWS_list_new.append(isPWS_list[i1])
                isFemale_list_new.append(isFemale_list[i1])
                SSI_reg.append(SSI4_list[i1])
            else:
                S_toDel.append(S_list[i1])
             
        if os.path.isdir(tbss_dir_new):
            print('Directory ' + tbss_dir_new + ' already exists.')
            bWrite = input('\tOverwrite? (0/1): ')
        else:
            bWrite = 1
    
        tbss_mat = tbss_mat.replace('.mat', '_corrSSI.mat')
        tbss_dir_0 = tbss_dir
        tbss_dir = tbss_dir_new
        S_list = S_list_new
        isPWS_list = isPWS_list_new
        isFemale_list = isFemale_list_new

        if bWrite == 1:
            print('Copying directory from ' + tbss_dir + ' to ' + tbss_dir_new)
            os.system('rm -rf ' + tbss_dir_new)
            os.system('cp -r ' + tbss_dir_0 + ' ' + tbss_dir_new)        
            os.system('rm -rf ' + os.path.join(tbss_dir_new, 'stats'))

            for std in S_toDel:
                print('Deleting the files belonging to mild PWS subject ' + std)
                os.system('rm -f ' + \
                          os.path.join(tbss_dir_new, 'origdata', std + '*'))
                os.system('rm -f ' + \
                          os.path.join(tbss_dir_new, 'FA', std + '*'))
        
            cwd = os.getcwd()
            os.chdir(tbss_dir_new)
            os.system('tbss_3_postreg -S')
            os.system('tbss_4_prestats %f' % SKEL_FA_THRESH)
            os.chdir(cwd)
    else:
        S_list_new = S_list
        isPWS_list_new = isPWS_list
        isFemale_list_new = isFemale_list
        S_toDel = []

    # === Write to design matrix file ===
    if not bCorrSSI:
        if bIncGen or bIncAge:
            mat_txt = '/NumWaves 3\n'
        else:
            mat_txt = '/NumWaves 2\n'

        mat_txt += '/NumPoints ' + str(len(S_list)) + '\n'
        mat_txt += '/PPheights 1 1\n'
        mat_txt += '/Matrix\n'
        for i1 in range(len(S_list)):
            if bIncGen:
                mat_txt += '1 ' + str(isPWS_list[i1]) + \
                           ' ' + str(isFemale_list[i1]) + '\n'
            elif bIncAge:
                mat_txt += '1 ' + str(isPWS_list[i1]) + \
                           ' ' + str(age_list[i1]) + '\n'
            else:
                mat_txt += '1 ' + str(isPWS_list[i1]) + '\n'
    else:    # Write the mat for corrSSI
        if bIncGen:
            mat_txt = '/NumWaves 3\n'
        else:
            mat_txt = '/NumWaves 2\n'

        SSI_mean_list = [-1 * sum(SSI_reg) / len(SSI_reg)] * len(SSI_reg)
        SSI_reg_dm = [sum(pair) for pair in zip(SSI_reg, SSI_mean_list)]
        mat_txt += '/NumPoints ' + str(len(S_list)) + '\n'
        mat_txt += '/PPheights 1 1\n'
        mat_txt += '/Matrix\n'
        for i1 in range(len(S_list)):
            if bIncGen:
                mat_txt += '1 ' + str(SSI_reg_dm[i1]) + \
                           ' ' + str(isFemale_list[i1]) + '\n'
            elif bIncAge:
                mat_txt += '1 ' + str(SSI_reg_dm[i1]) + \
                           ' ' + str(age_list[i1]) + '\n'
            else:
                mat_txt += '1 ' + str(SSI_reg_dm[i1]) + '\n'

    if os.path.isfile(tbss_mat):
        print('File ' + tbss_mat + ' already exists.')
        bWrite = input('\tOverwrite? (0/1): ')
    else:
        bWrite = 1

    if bWrite == 1:
        tbss_mat_f = open(tbss_mat, 'w')
        tbss_mat_f.write(mat_txt)
        tbss_mat_f.close()
        print('INFO: Written to ' + tbss_mat)


    # === Write to contrast vector file === 
    if bIncGen or bIncAge:
        con_txt = '/NumWaves 3\n'
        con_txt += '/NumContrasts 1\n'
        con_txt += '/PPheights 1 1 1\n'
    else:
        con_txt = '/NumWaves 2\n'
        con_txt += '/NumContrasts 1\n'
        con_txt += '/PPheights 1 1\n'

    con_txt += '/Matrix\n'
    con_txt_rev = con_txt

    if bIncGen or bIncAge:
        con_txt += '0 -1 0\n'
        con_txt_rev += '0 1 0\n'
    else:
        con_txt += '0 -1\n'
        con_txt_rev += '0 1\n'

    if os.path.isfile(tbss_con):
        print('File ' + tbss_con + ' already exists.')
        bWrite = input('\tOverwrite? (0/1): ')
    else:
        bWrite = 1

    if bWrite == 1:
        tbss_con_f = open(tbss_con, 'w')
        tbss_con_f.write(con_txt)
        tbss_con_f.close()
        print('INFO: Written to ' + tbss_con)

        tbss_con_rev_f = open(tbss_con_rev, 'w')
        tbss_con_rev_f.write(con_txt_rev)
        tbss_con_rev_f.close()
        print('INFO: Written to ' + tbss_con_rev)

    # === Copy directory structure and remove subjects if necessary ===
    """
    if bIncGen or bIncAge or bMaleOnly:
        if os.path.isdir(tbss_dir_new):
            print('File ' + tbss_dir_new + ' already exists.')
            bWrite = input('\tOverwrite? (0/1): ')
        else:
            bWrite = 1
    
        if bWrite == 1:
            print('Copying directory from ' + tbss_dir + ' to ' + tbss_dir_new)
            os.system('rm -rf ' + tbss_dir_new)
            os.system('cp -r ' + tbss_dir + ' ' + tbss_dir_new)
            os.system('rm -rf ' + os.path.join(tbss_dir_new, 'stats'))

            for std in S_toDel:
                print('Deleting the files belonging to mild PWS subject ' + std)
                os.system('rm ' + os.path.join(tbss_dir_new, 'origdata', std + '*'))
                os.system('rm ' + os.path.join(tbss_dir_new, 'FA', std + '*'))

            cwd = os.getcwd()
            os.chdir(tbss_dir_new)
            os.system('tbss_3_postreg -S')
            os.system('tbss_4_prestats 0.2')
            os.chdir(cwd)
        
            if bNoMild:
                tbss_mat = tbss_mat.replace('.mat', '_noMild.mat')
            tbss_dir = tbss_dir_new
            S_list = S_list_new
            isPWS_list = isPWS_list_new
            isFemale_list = isFemale_list_new
        else:
            sys.exit(1)

    sys.exit(0)
    """

    ## Offer to run the randomise
    skel = os.path.join(tbss_dir, 'stats', 'all_FA_skeletonised.nii.gz')
    mask = os.path.join(tbss_dir, 'stats', 'mean_FA_skeleton_mask.nii.gz')
    check_file(skel)
    check_file(mask)

    if bIncGen:
        if not bCorrSSI:
            output_prefix = \
                os.path.join(tbss_dir, 'stats', 'tbss_bgc_incGen')
        else:
            output_prefix = \
                os.path.join(tbss_dir, 'stats', 'tbss_corrSSI_incGen')
    elif bIncAge:
        if not bCorrSSI:
            output_prefix = \
                os.path.join(tbss_dir, 'stats', 'tbss_bgc_incAge')
        else:
            output_prefix = \
                os.path.join(tbss_dir, 'stats', 'tbss_corrSSI_incAge')
    else:
        if not bCorrSSI:
            output_prefix = os.path.join(tbss_dir, 'stats', 'tbss_bgc')
        else:
            output_prefix = os.path.join(tbss_dir, 'stats', 'tbss_corrSSI')

    randomise_500_cmd = 'randomise -i ' + skel + ' -o ' + \
                         output_prefix + '_500 ' + ' -m ' + mask + \
                         ' -d ' + tbss_mat + ' -t ' + tbss_con + \
                         ' -n 500 ' + '--T2 -V'
    randomise_1e4_cmd = 'randomise -i ' + skel + ' -o ' + \
                        output_prefix + '_1e4 ' + ' -m ' + mask + \
                        ' -d ' + tbss_mat + ' -t ' + tbss_con + \
                        ' -n 10000 ' + '--T2 -V'
    randomise_500_rev_cmd = 'randomise -i ' + skel + \
                            ' -o ' + output_prefix + '_rev_500 ' + \
                            ' -m ' + mask + ' -d ' + tbss_mat + \
                            ' -t ' + tbss_con_rev + ' -n 500 ' + '--T2 -V'
    randomise_1e4_rev_cmd = 'randomise -i ' + skel + \
                            ' -o ' + output_prefix + '_rev_1e4 ' + \
                            ' -m ' + mask + ' -d ' + tbss_mat + \
                            ' -t ' + tbss_con_rev + ' -n 10000 ' + '--T2 -V'

    print('\nrandomise_500_cmd = ' + randomise_500_cmd + '\n')
    print('\nrandomise_1e4_cmd = ' + randomise_1e4_cmd + '\n')
    print('\nrandomise_500_rev_cmd = ' + randomise_500_rev_cmd + '\n')
    print('\nrandomise_1e4_rev_cmd = ' + randomise_1e4_rev_cmd + '\n')

    jobName = 'rdms';
    if bIncGen:
        jobName += '_incGen'
    if bIncAge:
        jobName += '_incAge'
    if bNoMild:
        jobName += '_noMild'

    bRun = input('Run randomise commands? (0 - No / 1 - 500 / 2 - 500&10000): ')
    if bRun == 2:
        os.system('ezsub -q max30 -c "' + randomise_500_cmd + \
                  '" -n ' + jobName + '_500')
        os.system('ezsub -q max30 -c "' + randomise_1e4_cmd + \
                  '" -n ' + jobName + '_1e4')
        os.system('ezsub -q max30 -c "' + randomise_500_rev_cmd + \
                  '" -n ' + jobName + '_rev_500')
        os.system('ezsub -q max30 -c "' + randomise_1e4_rev_cmd + \
                  '" -n ' + jobName + '_rev_1e4')
    elif bRun == 1:
        os.system('ezsub -q max30 -c "' + randomise_500_cmd + '" -n ' + jobName + '_500')
        os.system('ezsub -q max30 -c "' + randomise_500_rev_cmd + '" -n ' + jobName + '_rev_500')
