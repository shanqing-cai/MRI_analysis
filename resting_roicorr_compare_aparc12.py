#!/usr/bin/python

import os
import sys
import glob
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat

from get_qdec_info import get_qdec_info
from fdr import fdr
from exclusion_paradigms import exclParad

qdec_fn = '/users/cais/STUT/FSDATA/qdec/qdec.table.dat'
bips_resting_dir = '/users/cais/STUT/analysis/'
signifThresh = 0.0001
#signifThresh = 0.01

seedROI = ''
corrMeas = ''
bFDR = False
FDR = np.nan


#sys.exit()

def get_sIDs(qdec_fn, bips_resting_dir, resfn):
    sIDs = {'PWS': [], 'PFS': []}

    qdec_f = open(qdec_fn, 'r')
    qdec_txt = qdec_f.read()
    qdec_f.close()
    qdec_txt = qdec_txt.split('\n')

    for t_line in qdec_txt:
        if len(t_line) == 0:
            continue

        t_lin = t_line.replace('\t', ' ')
        t_lin = t_lin.split(' ')
        if t_lin[0] == 'fsid':
            continue
        if len(t_lin) == 0:
            continue

        while t_lin.count('') > 0:
            t_lin.remove('')
        
        t_dgn = t_lin[2]
        t_sID = t_lin[0]

        hemi1 = 'lh'
        hemi2 = 'rh'
        if os.path.isfile(os.path.join(bips_resting_dir, \
                          t_sID, resfn)) and \
           os.path.isfile(os.path.join(bips_resting_dir, \
                          t_sID, resfn)):
            sIDs[t_dgn].append(t_sID)
        else:
            print('WARNING: %s files not found for subject %s'%(resfn, t_sID))
    
    return sIDs

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: resting_roicorr_compare_aparc12.py <mode> [options]')
        print('\tmode = {fullspectrum, z_no_outliers_bandpassed, znbp2, bpnrm}')
        sys.exit(0)

    if len(sys.argv) > 1: # Additional input arguments
        if sys.argv.count('--seedROI') == 1:
            seedROI = sys.argv[sys.argv.index('--seedROI') + 1]
            print('seedROI = %s\n'%seedROI)
            #sys.exit(0)
        if sys.argv.count('--corr') == 1:
            corrMeas = sys.argv[sys.argv.index('--corr') + 1]
            print('corrMeas = %s\n'%corrMeas)
        if sys.argv.count('--fdr') == 1:
            bFDR = True
            FDR = float(sys.argv[sys.argv.index('--fdr') + 1])
            print('bFDR = True')
            print('FDR = %f'%float(FDR))


    mode = sys.argv[1]
    if not (mode == 'fullspectrum' or mode == 'z_no_outliers_bandpassed' \
            or mode == "bpnrm" or mode == "znbp2"):
        raise ValueError, 'Unrecognized mode: %s'%mode
        sys.exit(0)

    if mode == 'fullspectrum':
        resfn = 'roi_corr_aparc12.csc.fullspectrum.pkl'
    elif mode == 'z_no_outliers_bandpassed':
        resfn = 'roi_corr_aparc12.csc.z_no_outliers_bandpassed.pkl'
    elif mode == 'znbp2':
        resfn = 'roi_corr_aparc12.csc.z_no_outliers_bandpassed2.pkl'
    elif mode == "bpnrm":
        resfn = 'roi_corr_aparc12.csc.bpnrm.pkl'

    sIDs = get_sIDs(qdec_fn, bips_resting_dir, resfn)

    #--- Apply exclusion paradigm ---#
    exclName = ''
    for arg in sys.argv:
        if arg.startswith('--exclParad:'):
            exclName = arg.replace('--exclParad:', '')    

    grps = sIDs.keys()
    if exclName != '':
        exclS = exclParad[exclName]
        
        for grp in grps:
            for e_sID in exclS:
                if sIDs[grp].count(e_sID) == 1:
                    sIDs[grp].remove(e_sID)
                    print('Removing subject %s from group %s'%(e_sID, grp))

        print('')
    #--- ~Apply exclusion paradigm ---#

    print("N(PWS) = %d"%(len(sIDs["PWS"])))
    print("N(PFS) = %d"%(len(sIDs["PFS"])))
    print("")

    res = np.load(os.path.join(bips_resting_dir, sIDs['PFS'][0], resfn))
    nROIs = len(res['b_rois'])
    b_rois = res['b_rois']
    b_ids = res['b_ids']    

    p_tab = np.array([[np.nan] * nROIs] *nROIs)

    nPFS = len(sIDs['PFS'])
    nPWS = len(sIDs['PWS'])
    
    corr_z_tab = {'PFS': np.array([[np.nan] * nROIs * nROIs] * nPFS), \
                'PWS': np.array([[np.nan] * nROIs * nROIs] * nPWS)}
    SSI4 = np.array([np.nan] * nPWS)
    meas = {'PFS': np.array([np.nan] * nPFS), \
            'PWS': np.array([np.nan] * nPWS)}

    for grp in corr_z_tab.keys():
        print('Processing data from group %s...'%grp)
        for (i1, sID) in enumerate(sIDs[grp]):
            print('Processing data from subject %s...'%sID)

            res = np.load(os.path.join(bips_resting_dir, sID, resfn))
            t_zs = np.matrix.reshape(res['corr_tab'], [1, nROIs * nROIs])[0]
            
            # Fisher (z) transformation
            corr_z_tab[grp][i1] = np.log((1 + t_zs) / (1 - t_zs)) * 0.5
            # sys.exit(0)

            if grp == 'PWS':
                t_SSI = get_qdec_info(sID, 'SSI')
                SSI4[i1] = float(t_SSI)

            if len(corrMeas) > 0:
                t_meas = get_qdec_info(sID, corrMeas)
                try:
                    meas[grp][i1] = float(t_meas)
                except:
                    print('\tWARNING: %s measure for %s is not valid'%(corrMeas, sID))
                       
        print('')
        corr_z_tab[grp] = np.matrix.transpose(corr_z_tab[grp])

    ### Write to a clean mat file ###
    idx_keep = []
    for i0 in range(len(corr_z_tab['PFS'])):
        if ~np.isnan(corr_z_tab['PFS'][i0][0]):
            idx_keep.append(i0)

    corr_z = {'PFS': [], 'PWS': []}
    corr_z['PFS'] = corr_z_tab['PFS'][idx_keep]
    corr_z['PWS'] = corr_z_tab['PWS'][idx_keep]

    a_roi1 = []
    a_roi2 = []
    kidx = 0
    for roi1 in res['b_rois']:
        for roi2 in res['b_rois']:
            if idx_keep.count(kidx) == 1:
                a_roi1.append(roi1)
                a_roi2.append(roi2)
            kidx = kidx + 1

    if exclName == '':
        saveMatFN = os.path.join(bips_resting_dir, 'corr_z_roi_aparc12.%s.mat'%mode)
    else:
        saveMatFN = os.path.join(bips_resting_dir, \
                                 'corr_z_roi_aparc12.%s.%s.mat'%(mode, exclName))

    savemat(saveMatFN, \
            {'corr_z': corr_z, 'a_roi1': a_roi1, 'a_roi2': a_roi2, \
             'sIDs': sIDs})
    print('\nSaved two-group matrix to mat file %s\n'%saveMatFN)
        
    ### ~Write to a clean mat file ###

    #sys.exit()
    meas_all = np.concatenate((meas['PFS'], meas['PWS']))
    
    showZScores = {'PFS': [], 'PWS': []}
    sigROIs = {'roi1': [], 'roi2': []}
    t_showZScores = {'PFS': [], 'PWS': []}
    t_sigROIs = {'roi1': [], 'roi2': []}

    showZScores_PWS_wSSI4 = []
    t_showZScores_PWS_wSSI4 = []

    showZScores_opt = []
    showMeas_opt = []

    k = 0
    nTotComps = 0
    nTotCorrs_wSSI4 = 0
    nTotCorrs_opt = 0
    nSigDiff = 0
    all_ps = np.array([])
    all_ps_wSSI4 = np.array([])
    
    #r_corr_w_SSI4 = []
    #p_corr_w_SSI4 = []
    nSigCorrs_wSSI4 = 0
    sigROIs_wSSI4 = {'roi1': [], 'roi2': []}
    t_sigROIs_wSSI4 = {'roi1': [], 'roi2': []}

    rSigCorr_wSSI4 = []
    t_rSigCorr_wSSI4 = []
    pSigCorr_wSSI4 = []
    t_pSigCorr_wSSI4 = []

    nSigCorrs_opt = 0
    sigROIs_opt = {'roi1': [], 'roi2': []}
    rSigCorr_opt = []
    pSigCorr_opt = []

    for i1 in range(nROIs):
        roi1 = b_rois[i1]

        for i2 in range(nROIs):
            roi2 = b_rois[i2]

            if len(seedROI) > 0:
                if not ((roi1 == seedROI) or (roi2 == seedROI)):
                    k = k + 1
                    continue

            z_PFS = corr_z_tab['PFS'][k]
            z_PWS = corr_z_tab['PWS'][k]

            #t_r_corr_w_SSI4 = np.corrcoef(z_PWS, SSI4)
            #t_r_corr_w_SSI4 = t_r_corr_w_SSI4[0][1]
            #sys.exit(0)
            
            if len(np.where(np.isnan(z_PFS))[0]) > 0:
                k = k + 1
                continue;

            nTotCorrs_wSSI4 = nTotCorrs_wSSI4 + 1
            (t_r, t_p) = stats.pearsonr(z_PWS, SSI4)
            #if t_p < signifThresh:
            #nSigCorrs_wSSI4 = nSigCorrs_wSSI4 + 1
            t_sigROIs_wSSI4['roi1'].append(roi1)
            t_sigROIs_wSSI4['roi2'].append(roi2)
            t_rSigCorr_wSSI4.append(t_r)
            t_pSigCorr_wSSI4.append(t_p)
            t_showZScores_PWS_wSSI4.append(z_PWS)

            if len(corrMeas) > 0:
                nTotCorrs_opt = nTotCorrs_opt + 1
                z_all = np.concatenate((z_PFS, z_PWS))
                idx_keep = ~np.isnan(meas_all)
                (t_r, t_p) = stats.pearsonr(z_all[idx_keep], meas_all[idx_keep])
               
                if t_p < signifThresh:
                    nSigCorrs_opt = nSigCorrs_opt + 1
                    sigROIs_opt['roi1'].append(roi1)
                    sigROIs_opt['roi2'].append(roi2)
                    rSigCorr_opt.append(t_r)
                    pSigCorr_opt.append(t_p)
                    showZScores_opt.append(z_all[idx_keep])
                    showMeas_opt.append(meas_all[idx_keep])
            

            nTotComps = nTotComps + 1
            (t, p) = stats.ttest_ind(z_PFS, z_PWS)
            p_tab[i1][i2] = p
            all_ps = np.append(all_ps, p)
            
            #if p < signifThresh:
            #    nSigDiff = nSigDiff + 1
            if np.mean(z_PWS) > np.mean(z_PFS):
                strDirection = 'PWS > PFS'
            else:
                strDirection = 'PWS < PFS'
            #print('ROI1 = %s, ROI2 = %s: p = %f (%s)'%(b_rois[i1], b_rois[i2], p, strDirection))
            #print('ROI1 = %s, ROI2 = %s: p = %f (%s)'%(roi1, roi2, p, strDirection))
            t_showZScores['PFS'].append(z_PFS)
            t_showZScores['PWS'].append(z_PWS)
            t_sigROIs['roi1'].append(roi1)
            t_sigROIs['roi2'].append(roi2)

            k = k + 1

    
    

    print('Between-group comparisons of correlation z-scores')
    if bFDR:
        signifThresh = fdr(all_ps, FDR)
        print('**FDR** (q = %f) --> signifThresh = %f\n'%(FDR, signifThresh))
    else:
        print('Pre-specified signifThresh = %f\n'%signifThresh);

    ### Determine the significant differences ###
    keep_idx = []
    for i1 in range(nTotComps):
        if all_ps[i1] < signifThresh:
            keep_idx.append(i1)
    nSigDiff = len(keep_idx)

    for idx in keep_idx:
        p = all_ps[idx]
        zScores_PFS = t_showZScores['PFS'][idx]
        zScores_PWS = t_showZScores['PWS'][idx]
        roi1 = t_sigROIs['roi1'][idx]
        roi2 = t_sigROIs['roi2'][idx]

        if np.mean(zScores_PWS) > np.mean(zScores_PFS):
            strDirection = 'PWS > PFS'
        else:
            strDirection = 'PWS < PFS'
        print('ROI1 = %s, ROI2 = %s: p = %f (%s)'%(roi1, roi2, p, strDirection))


        showZScores['PFS'].append(zScores_PFS)
        showZScores['PWS'].append(zScores_PWS)
        sigROIs['roi1'].append(roi1)
        sigROIs['roi2'].append(roi2)

    ### Visualize the significant between-group differences in z ###
    print('Total number of comparisons: %d'%nTotComps)
    print('Number of significant differences: %d'%nSigDiff)

    plts = []
    xlim_left = 0.5
    xlim_right = 2.5
    ylim_bottom = -0.5
    ylim_top = 2.
    for i1 in range(nSigDiff):
        plts.append(plt)
        n_PFS = len(showZScores['PFS'][i1])
        n_PWS = len(showZScores['PWS'][i1])
        mean_PFS = np.mean(showZScores['PFS'][i1])
        mean_PWS = np.mean(showZScores['PWS'][i1])
        ste_PFS = np.std(showZScores['PFS'][i1]) / n_PFS
        ste_PWS = np.std(showZScores['PFS'][i1]) / n_PWS
        
        #plts[-1].ion()
        plts[-1].figure()
        plts[-1].plot([xlim_left, xlim_right], [0, 0])
        plts[-1].plot([1] * n_PFS, showZScores['PFS'][i1], 'ko')
        plts[-1].plot([2] * n_PWS, showZScores['PWS'][i1], 'ro')
        plts[-1].plot(1.2, mean_PFS)
        plts[-1].plot(2.2, mean_PWS)
        plts[-1].plot([1.2, 1.2], [mean_PFS - ste_PFS, mean_PFS + ste_PFS], 'k-')
        plts[-1].plot([2.2, 2.2], [mean_PWS - ste_PWS, mean_PWS + ste_PWS], 'r-')
        plts[-1].xlim([xlim_left, xlim_right])
        plts[-1].ylim([ylim_bottom, ylim_top])
        plts[-1].title(sigROIs['roi1'][i1] + ' x ' + sigROIs['roi2'][i1])
        plts[-1].xticks([1, 2])
        plts[-1].ylabel('Fisher transformed corr. coeff. (z)')
        plts[-1].show()
    
    ### Determine significant correlations with SSI4 ###
    print('\n\nCorrelation with SSI4 scores in the PWS group')
    all_ps_wSSI4 = np.array(t_pSigCorr_wSSI4)
    if bFDR:
        signifThresh_wSSI4 = fdr(all_ps_wSSI4, FDR)
        print('**FDR** (q = %f) --> signifThresh = %f\n'%(FDR, signifThresh_wSSI4))
    else:
        signifThresh_wSSI4 = signifThresh
        print('Pre-specified signifThresh = %f\n'%signifThresh_wSSI4);

    keep_idx = []
    for i1 in range(len(all_ps_wSSI4)):
        if all_ps_wSSI4[i1] < signifThresh_wSSI4:
            keep_idx.append(i1)
    nSigCorrs_wSSI4 = len(keep_idx)

    for idx in keep_idx:
        p = all_ps_wSSI4[idx]
        r = t_rSigCorr_wSSI4[idx]
        zScores = t_showZScores_PWS_wSSI4[idx]
        roi1 = t_sigROIs_wSSI4['roi1'][idx]
        roi2 = t_sigROIs_wSSI4['roi2'][idx]

        print('ROI1 = %s, ROI2 = %s: p = %f, r = %f'%(roi1, roi2, p, r))

        showZScores_PWS_wSSI4.append(zScores)
        sigROIs_wSSI4['roi1'].append(roi1)
        sigROIs_wSSI4['roi2'].append(roi2)

    print('Total number of correlations with SSI4: %d'%nTotCorrs_wSSI4)
    print('Number of significant correlations with SSI4 (PWS): %d'%nSigCorrs_wSSI4)

    plts_SSI4 = []
    for i1 in range(nSigCorrs_wSSI4):
        plts_SSI4.append(plt)
        n_PWS = len(showZScores_PWS_wSSI4[i1])

        plts_SSI4[-1].figure()
        plts_SSI4[-1].plot(SSI4, showZScores_PWS_wSSI4[i1], 'ro')
        plts_SSI4[-1].xlabel('SSI4 score')
        plts_SSI4[-1].ylabel('z score (%s x %s)'%(sigROIs_wSSI4['roi1'][i1], sigROIs_wSSI4['roi2'][i1]))
        plts_SSI4[-1].show()

        #print('r = %f, p = %f @ %s x %s' \
        #          %(rSigCorr_wSSI4[i1], pSigCorr_wSSI4[i1], sigROIs_wSSI4['roi1'][i1], sigROIs_wSSI4['roi2'][i1]))
    

    # Visualize the significant correlations (with SSI4)
    #''' 
    

    
    #'''
    
    # Visualize the significant correlations (optional)
    print('\nTotal number of correlations with %s: %d'%(corrMeas, nTotCorrs_opt))
    print('Number of significant correlations with %s (all Ss): %d'%(corrMeas, nSigCorrs_opt))
    
    if len(corrMeas) > 0:
        plts_opt = []
        for i1 in range(nSigCorrs_opt):
            plts_opt.append(plt)
            
            plts_opt[-1].figure()
            plts_opt[-1].plot(showMeas_opt[i1], showZScores_opt[i1], 'bo')
            plts_opt[-1].xlabel('Measure: %s'%corrMeas)
            plts_opt[-1].ylabel('z score (%s x %s)'%(sigROIs_opt['roi1'][i1], sigROIs_opt['roi2'][i1]))
            plts_opt[-1].show()

            print('r = %f, p = %f @ %s x %s' \
                  %(rSigCorr_opt[i1], pSigCorr_opt[i1], sigROIs_opt['roi1'][i1], sigROIs_opt['roi2'][i1]))
            
