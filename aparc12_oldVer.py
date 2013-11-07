#!/usr/bin/python
DEBUG = False

import numpy as np

# 3rd column: N - non-speech; S - speech
aROIs = [   ['FP',     'Prefrontal',   'N'], \
            ['aMFg',   'Prefrontal',   'N'], \
            ['aIFs',   'Prefrontal',   'N'], \
            ['SFg',    'Prefrontal',   'N'], \
            ['aFO',    'Prefrontal',   'S'], \
            ['FOC',    'Prefrontal',   'N'], \
            ['pMFg',   'Prefrontal',   'N'], \
            ['pIFs',   'Prefrontal',   'S'], \
            ['dIFt',   'Prefrontal',   'N'], \
            ['dIFo',   'Prefrontal',   'S'], \
            ['vIFt',   'Prefrontal',   'N'], \
            ['vIFo',   'Prefrontal',   'S'], \
            ['pFO',    'Prefrontal',   'S'], \
            ['FMC',    'Prefrontal',   'N'], \
            ['aINS',   'Insular',      'S'], \
            ['pINS',   'Insular',      'S'], \
            ['adPMC',  'Premotor',     'N'], \
            ['mdPMC',  'Premotor',     'N'], \
            ['pdPMC',  'Premotor',     'S'], \
            ['midPMC', 'Premotor',     'S'], \
            ['vPMC',   'Premotor',     'S'], \
            ['SMA',    'Premotor',     'S'], \
            ['preSMA', 'Premotor',     'S'], \
            ['aCO',    'Precentral',   'S'], \
            ['dMC',    'Precentral',   'S'], \
            ['midMC',  'Precentral',   'S'], \
            ['vMC',    'Precentral',   'S'], \
            ['pCO',    'Postcentral',  'S'], \
            ['dSC',    'Postcentral',  'S'], \
            ['vSC',    'Postcentral',  'S'], \
            ['aSMg',   'PPC',          'S'], \
            ['pSMg',   'PPC',          'N'], \
            ['PO',     'PPC',          'S'], \
            ['SPL',    'PPC',          'N'], \
            ['Ag',     'PPC',          'N'], \
            ['PCN',    'PPC',          'N'], \
            ['TP',     'Temporal',     'N'], \
            ['PP',     'Temporal',     'S'], \
            ['H',     'Temporal',     'S'], \
            ['PT',     'Temporal',     'S'], \
            ['aSTg',   'Temporal',     'S'], \
            ['pSTg',   'Temporal',     'S'], \
            ['pdSTs',  'Temporal',     'S'], \
            ['adSTs',  'Temporal',     'N'], \
            ['pvSTs',  'Temporal',     'N'], \
            ['avSTs',  'Temporal',     'N'], \
            ['aMTg',   'Temporal',     'N'], \
            ['pMTg',   'Temporal',     'N'], \
            ['pITg',   'Temporal',     'N'], \
            ['aCG',    'Cingulate',    'S'], \
            ['midCG',  'Cingulate',    'N'], \
            ['pCG',    'Cingulate',    'N'], \
            ['OC',     'Occipital',    'N'], \
            ['MTO',    'Occipital',    'N'], \
            ['ITO',    'Occipital',    'N'], \
            ['Lg',     'Occipital',    'N'], \
            ['pPH',    'Parahippocampal',  'N'], \
            ['aPH',    'Parahippocampal',  'N'], \
            ['SCC',    'Parahippocampal',  'N']]


aROIs = np.array(aROIs)
nROIs = len(aROIs)

if DEBUG:
    print("len(aROIs) = %d" % len(aROIs)) # DEBUG

roiNames = aROIs[:, 0]

if DEBUG:
    print("len(np.unique(roiNames)) = %d" % len(np.unique(roiNames))) # DEBUG

lobeNames = aROIs[:, 1]
uLobeNames = np.unique(lobeNames)

if DEBUG:
    print("len(np.unique(lobeNames) = %d" % len(uLobeNames)) # DEBUG

isSpeech = aROIs[:, 2]
nSpeechROIs = len(np.nonzero(isSpeech == 'S')[0])

if DEBUG:
    print("nSpeechROIs = %d" % nSpeechROIs) # DEBUG

activROI_uc_fn_wc = "/users/cais/STUT/scripts/activROIs_uc_thr%.3f.mat"

from scai_utils import *
from scipy.io import loadmat

def get_aparc12_cort_rois(lobe="all", bSpeech=False):
#if __name__ == "__main__":
#    lobe = "all"
#    bSpeech = "speech_2g_lh_0.01"

    if lobe == "all":
        idxLobe = np.array(range(nROIs))
    else:
        if (len(np.nonzero(uLobeNames == lobe)[0])) == 0:
            idxLobe = np.array([])
            raise Exception, "Unrecognized lobe name: %s" % lobe
        else:
            idxLobe = np.nonzero(lobeNames == lobe)[0]
            # ROIs = roiNames[idx]

    if bSpeech == True:
        idxSpeech = np.nonzero(isSpeech == "S")[0]
    elif bSpeech == False:
        idxSpeech = np.array(range(nROIs))
    elif bSpeech.startswith("speech_PFS_lh") \
            or bSpeech.startswith("speech_PFS_rh") \
            or bSpeech.startswith("speech_PWS_lh") \
            or bSpeech.startswith("speech_PWS_rh") \
            or bSpeech.startswith("speech_2g_lh") \
            or bSpeech.startswith("speech_2g_rh"):
        
        if bSpeech.count('_') != 3:
            raise Exception, "Number of underlines does not equal 3"

        t_grp = bSpeech.split('_')[1]
        t_hemi = bSpeech.split('_')[2]
        t_thr = float(bSpeech.split('_')[3])

        activROI_uc_fn = activROI_uc_fn_wc % t_thr
        check_file(activROI_uc_fn)

        roiSet = loadmat(activROI_uc_fn)
        roiSet = roiSet['roiSet']
        
        # print(t_grp)
        # print(t_hemi)

        hemis = ['lh', 'rh']
        hemii = hemis.index(t_hemi)
        if t_grp == "PFS" or t_grp == "PWS":
            rois = roiSet[t_grp][0][0][0][0][hemii][0]
        else:
            rois_PFS = roiSet["PFS"][0][0][0][0][hemii][0]
            rois_PWS = roiSet["PWS"][0][0][0][0][hemii][0]
            rois = np.union1d(rois_PFS, rois_PWS)

        rois = list(rois)
        for (i0, t_roi) in enumerate(rois):
            t_roi = str(t_roi)[3 : -2]
            rois[i0] = t_roi

        a_rois = list(aROIs[:, 0])
        idxSpeech = []
        for (i0, t_roi) in enumerate(a_rois):
            if rois.count(t_roi) == 1:
                idxSpeech.append(i0)

    else:
        raise Exception, "Unrecognized mode: %s" % bSpeech
                
    idx = np.intersect1d(idxLobe, idxSpeech)
    
    if DEBUG:
        print(idxLobe)
        print(idxSpeech)
        print(idx)
    ROIs = list(roiNames[idx])

    return ROIs
