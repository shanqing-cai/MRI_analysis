#!/usr/bin/python

import os
import sys
import xml.etree.ElementTree as ET

def get_bad_frames(xmlfn, nFrames=70):
    if not os.path.isfile(xmlfn):
        raise Exception, "%s: file doesn't exist: %s" % \
                         (get_bad_grad_dirs.__name__, xmlfn)

    tree = ET.parse(xmlfn)
    root = tree.getroot()

    grad_cnt = 0
    qc_indices = []

    for child0 in root:
        #print child0.attrib["parameter"]        
        if child0.attrib["parameter"] == "DWI Check":
            for (i1, child1) in enumerate(child0):
                #print child1.attrib
                if child1.attrib["parameter"].startswith("gradient_"):
                    grad_cnt += 1
                    for (i2, child2) in enumerate(child1):
                        if len(child2.attrib) == 0:
                            continue
                        if child2.attrib["parameter"] == "QC_Index":
                            qc_indices.append(int(child2[0].text))

    bad_frame_indices = []
    for (i0, t_idx) in enumerate(qc_indices):
        if t_idx < 0:
            bad_frame_indices.append(i0)
    

    if grad_cnt != nFrames:
        raise Exception, "%s, Erroneous number of DWI Check entries in xml file: %s" % (function.__name__, xmlfn)
    
        
    #print("grad_cnt = %d\n" % grad_cnt)
    return bad_frame_indices

if __name__ == "__main__":
    xmlfn = "/users/cais/STUT/analysis/dti/S39/S39_XMLQCResult.xml"
    
    bad_frame_idx = get_bad_frames(xmlfn)

    
