#!/usr/bin/python

import sys
import os
import glob

## CONFIG: Paths
qdec_dat = '/users/cais/STUT/FSDATA/qdec/qdec.table.dat'
qdec_comp = '/users/cais/STUT/FSDATA/qdec/qdec_PWS.dat'

## Get the list of subjects
if len(sys.argv) == 1:
    sys.exit('Wrong syntax')

tbss_dir = sys.argv[1]

bMaleOnly = False
bIncGen = False
bNoMild = False
bCorrSSI = False

if len(sys.argv) > 2:
    for opt in sys.argv[2:]:
        if opt == '-incGen':
            bIncGen = True
        elif opt == '-noMild':
            bNoMild = True
        elif opt == '-corrSSI':
            bCorrSSI = True
        elif opt == '-maleOnly':
            bMaleOnly = True
        else:
            sys.exit('Unrecognized option: ' + opt)

if bCorrSSI and bNoMild:
    sys.exit('Do NOT use the options -noMild and -corrSSI simultaneously.')

if bIncGen and bMaleOnly:
    sys.exit('Do NOT use the options -incGen and -maleOnly simultaneously.')

#if len(sys.argv) == 2 or len(sys.argv) == 3:
#    tbss_dir = sys.argv[1]
#    if tbss_dir[-1] == '/':
#        tbss_dir = tbss_dir[:-1]
#    if len(sys.argv) == 3:
#        incGenOpt = sys.argv[2]
#        if incGenOpt == '-incGen':
#            bIncGen = True
#        else:
#            sys.exit('Wrong syntax')

#if len(sys.argv) == 4:
#    noMild = sys.argv[2]
#    thresh_SSI4 = float(sys.argv[3])
#else:
thresh_SSI4 = 16

#if noMild == '-noMild':
#    bNoMild = True
#else:
#    bNoMild = False

if bIncGen:
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

if os.path.isdir(tbss_dir):
    print('tbss_dir = ' + tbss_dir)
else:
    sys.exit('ERROR: tbss_dir ' + tbss_dir + ' not found.')

orig_dir = os.path.join(tbss_dir, 'origdata')
if not os.path.isdir(orig_dir):
    sys.exit('ERROR: origdata directory not found')

if not os.path.isfile(qdec_comp):
    sys.exit('ERROR: qdec_comp file ' + qdec_comp + ' not found')

d_orig = glob.glob(os.path.join(orig_dir, 'S*.nii.gz'))
S_list = []
for fn in d_orig:
    S_list.append(fn.split('/')[-1].split('.')[0])
S_list.sort()
print('\nFound ' + str(len(S_list)) + ' subjects in ' + orig_dir + '\n')

isPWS_list = [0] * len(S_list)
isFemale_list = [0] * len(S_list)
done_list = [0] * len(S_list)
SSI4_list = [-1.0] * len(S_list)

##
if not os.path.isfile(qdec_dat):
    sys.exit('ERROR: qdec file: ' + qdec_dat + ' not found.')

qdec_f = open(qdec_dat, 'r')
t = qdec_f.read()
t = t.split('\n')
t = t[1 : ]

qdec_f.close()

for line in t:
    if len(t) == 0:
        continue

    S = line.split(' ')[0]

    try:
        idx = S_list.index(S)
    except:
        continue
    print("S = %s; idx = %d"%(S, idx))

    iPWS = -1
    try:
        iPWS = line.index('PWS')
    except:
        pass
    if iPWS != -1:
        isPWS_list[idx] = 1

    iFemale = -1
    try:
        iFemale = line.index('Female')
    except:
        pass
    if iFemale != -1:
        isFemale_list[idx] = 1

    done_list[idx] = 1

print('sum(done_list) = ' + str(sum(done_list)))
if sum(done_list) == len(S_list):
    print('\tOK')
else:
    sys.exit('ERROR: vector could not be generated for some of the subjects in the TBSS directory')

#sys.exit(0)

## Figure out the SSI4 scores for the PWS
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
            if t_sID == S_list[i1] or t_sID.replace('_500', '') == S_list[i1]:
                bFound = True
                SSI4_list[i1] = float(t_line_1[3])
                #sys.exit(0)
                break

        if not bFound:
            sys.exit('ERROR: SSI4 score not found for subject ' + S_list[i1])
        else:
            print(S_list[i1] + ': SSI4 = ' + str(SSI4_list[i1]))

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
    
    # Copy
    tbss_dir_new = tbss_dir + '_noMild'

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

    # Copy
    tbss_dir_new = tbss_dir + '_maleOnly'

elif bCorrSSI:      # Include SSI4 scores as a covariate
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
             
    # Copy
    tbss_dir_new = tbss_dir + '_corrSSI'
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
        os.system('tbss_4_prestats 0.2')
        os.chdir(cwd)


# Copy directory structure and remove subjects if necessary
if bIncGen or bMaleOnly:
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

## Write to file
if not bCorrSSI:
    if bIncGen:
        mat_txt = '/NumWaves 3\n'
    else:
        mat_txt = '/NumWaves 2\n'

    mat_txt += '/NumPoints ' + str(len(S_list)) + '\n'
    mat_txt += '/PPheights 1 1\n'
    mat_txt += '/Matrix\n'
    for i1 in range(len(S_list)):
        if bIncGen:
            mat_txt += '1 ' + str(isPWS_list[i1]) + ' ' + str(isFemale_list[i1]) + '\n'
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
            mat_txt += '1 ' + str(SSI_reg_dm[i1]) + ' ' + str(isFemale_list[i1]) + '\n'
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

#sys.exit(0)

if bIncGen:
    con_txt = '/NumWaves 3\n'
else:
    con_txt = '/NumWaves 2\n'
con_txt += '/NumContrasts 1\n'
con_txt += '/PPheights 1 1\n'
con_txt += '/Matrix\n'
con_txt_rev = con_txt

if bIncGen:
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

#sys.exit(0)

## Offer to run the randomise
skel = os.path.join(tbss_dir, 'stats', 'all_FA_skeletonised.nii.gz')
mask = os.path.join(tbss_dir, 'stats', 'mean_FA_skeleton_mask.nii.gz')

if os.path.isfile(skel) and os.path.isfile(mask):
    if bIncGen:
        if not bCorrSSI:
            output_prefix = os.path.join(tbss_dir, 'stats', 'tbss_bgc_incGen')
        else:
            output_prefix = os.path.join(tbss_dir, 'stats', 'tbss_corrSSI_incGen')
    else:
        if not bCorrSSI:
            output_prefix = os.path.join(tbss_dir, 'stats', 'tbss_bgc')
        else:
            output_prefix = os.path.join(tbss_dir, 'stats', 'tbss_corrSSI')
    randomise_500_cmd = 'randomise -i ' + skel + ' -o ' + output_prefix + '_500 ' + ' -m ' + mask + \
                                    ' -d ' + tbss_mat + ' -t ' + tbss_con + ' -n 500 ' + \
                                    '--T2 -V'
    randomise_1e4_cmd = 'randomise -i ' + skel + ' -o ' + output_prefix + '_1e4 ' + ' -m ' + mask + \
                                    ' -d ' + tbss_mat + ' -t ' + tbss_con + ' -n 10000 ' + \
                                    '--T2 -V'
    randomise_500_rev_cmd = 'randomise -i ' + skel + ' -o ' + output_prefix + '_rev_500 ' + ' -m ' + mask + \
                                    ' -d ' + tbss_mat + ' -t ' + tbss_con_rev + ' -n 500 ' + \
                                    '--T2 -V'
    randomise_1e4_rev_cmd = 'randomise -i ' + skel + ' -o ' + output_prefix + '_rev_1e4 ' + ' -m ' + mask + \
                                    ' -d ' + tbss_mat + ' -t ' + tbss_con_rev + ' -n 10000 ' + \
                                    '--T2 -V'

print('\nrandomise_500_cmd = ' + randomise_500_cmd + '\n')
print('\nrandomise_1e4_cmd = ' + randomise_1e4_cmd + '\n')
print('\nrandomise_500_rev_cmd = ' + randomise_500_rev_cmd + '\n')
print('\nrandomise_1e4_rev_cmd = ' + randomise_1e4_rev_cmd + '\n')

jobName = 'rdms';
if bIncGen:
    jobName += '_incGen'
if bNoMild:
    jobName += '_noMild'


bRun = input('Run randomise commands? (0/1/2): ')
if bRun == 2:
    os.system('ezsub -q max30 -c "' + randomise_500_cmd + '" -n ' + jobName + '_500')
    os.system('ezsub -q max30 -c "' + randomise_1e4_cmd + '" -n ' + jobName + '_1e4')
    os.system('ezsub -q max30 -c "' + randomise_500_rev_cmd + '" -n ' + jobName + '_rev_500')
    os.system('ezsub -q max30 -c "' + randomise_1e4_rev_cmd + '" -n ' + jobName + '_rev_1e4')
elif bRun == 1:
    os.system('ezsub -q max30 -c "' + randomise_500_cmd + '" -n ' + jobName + '_500')
    os.system('ezsub -q max30 -c "' + randomise_500_rev_cmd + '" -n ' + jobName + '_rev_500')
