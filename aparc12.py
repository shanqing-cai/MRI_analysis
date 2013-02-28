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
            ['AG',     'PPC',          'N'], \
            ['PCN',    'PPC',          'N'], \
            ['TP',     'Temporal',     'N'], \
            ['PP',     'Temporal',     'S'], \
            ['H',      'Temporal',     'S'], \
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
            ['LG',     'Occipital',    'N'], \
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


def get_aparc12_cort_rois(lobe="all", bSpeech=False):
    if lobe == "all":
        idxLobe = np.array(range(nROIs))
    else:
        if (len(np.nonzero(uLobeNames == lobe)[0])) == 0:
            idxLobe = np.array([])
            raise Exception, "Unrecognized lobe name: %s" % lobe
        else:
            idxLobe = np.nonzero(lobeNames == lobe)[0]
            # ROIs = roiNames[idx]

    if bSpeech:
        idxSpeech = np.nonzero(isSpeech == "S")[0]
    else:
        idxSpeech = np.array(range(nROIs))

    idx = np.intersect1d(idxLobe, idxSpeech)
    
    if DEBUG:
        print(idxLobe)
        print(idxSpeech)
        print(idx)
    ROIs = roiNames[idx]

    return ROIs

#    if bSpeech:
    """
    ROIs = ['H', 'PO', 'PP', 'PT', 'SMA', \
                'aCG', 'aCO', 'aFO', 'aINS', 'aSMg', \
                'dIFo', 'dMC', 'dSC', 'midMC', 'midPMC', \
                'pCO', 'pFO', 'pIFs', 'pINS', 'pSTg', \
                'pdPMC', 'pdSTs', 'preSMA', 'vIFo', 'vMC', \
                'vPMC', 'vSC']
    """
#    else:
        
        

    """
    ROIs = ['adPMC', 'pMFg', 'pIFs', 'dIFt', 'dIFo', 'vIFo', 'vIFt', 'pFO', 'aINS', 'mdPMC', \
        'pdPMC', 'midPMC', 'vPMC', 'aCO', 'pINS', 'dMC', 'midMC', 'vMC', \
        'pCO', 'dSC', 'vSC', 'aSMg', 'PO', \
        'TP', 'PP', 'H', 'PT', 'pSTg', 'pdSTs', 'adSTs', 'pvSTs', 'avSTs', \
        'aMTg', 'pMTg', \
        'SMA', 'preSMA', 'aCG', \
        'FP', 'aMFg', 'aIFs', 'SFg', 'aFO', 'FOC', 'SPL', 'pSMg', \
        'AG', 'OC', 'MTO', 'ITO', 'pITg', 'PCN', 'LG', 'pPH', \
        'aPH', 'midCG', 'pCG', 'SCC', 'FMC']
    """

