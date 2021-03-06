#!/usr/bin/env python

import os
import sys
import numpy as np
import nibabel as nb
import argparse
import mayavi.mlab as mlab

from scai_utils import *

COORD_FILE = "/home/cais/STUT/FSDATA/fsaverage2/mri/aparc12_roi_coords.txt"
STRUCT_VOL = "/home/cais/STUT/FSDATA/fsaverage2/mri/brain.nii.gz"
# STRUCT_VOL = "/home/cais/STUT/FSDATA/MNI152_1mm/mri/brain.nii.gz"

# --- To convert into input --- #
# componentFile = "corrSSI4_sigComponents_1.txt"
# edges = [[0, 1], [1, 2], [2, 3], [3, 0], [1, 3]]

DEFAULT_OPACITY=1.0

COMPONENT_CLRS = [(1.0, 0.75, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]

def translate_roi_name(roiName):
    if roiName.count("_Hg") > 0:
        return roiName.replace("_Hg", "_H")
    elif roiName.count("_aCGg") > 0:
        return roiName.replace("_aCGg", "_aCG")
    else:
        return roiName

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Render 3D network component image based on an input component file (from analyze_pt2_seedOnly_cmat.m. The thickness of the tubes are proportional to the sig value in the component file")
    ap.add_argument("componentFile", type=str)
    ap.add_argument("outputImgFN", type=str)
    ap.add_argument("--hemi", type=str, default="lh")
    ap.add_argument("--opacity", dest="t_opacity", \
                    type=float, default=DEFAULT_OPACITY, \
                    help="Opacity of 3D tubes (default=%f)" % DEFAULT_OPACITY)
    ap.add_argument("--noText", dest="bNoText", action="store_true", \
                    help="Do not plot the 3D text")
    ap.add_argument("--struct-vol", dest="altStructVol", type=str, 
                    help="Specify STRUCT_VOL that differs from the default: %s"% STRUCT_VOL)
    ap.add_argument("--coord-file", dest="altCoordFile", type=str, 
                    help="Specify coordinates file that differs from the default: %s" % COORD_FILE)
    ap.add_argument("--translate-roi-name", dest="bTranslateROIName", 
                    action="store_true", 
                    help="Translate ROI names for old-version of aparc12")
    ap.add_argument("--no-edges", dest="bNoEdges", 
                    help="Do not render the edges of the network")
    ap.add_argument("--vmax", dest="vmax", type=float, default=1600, 
                    help="vmax for brain volume rendering (default: 1600)")
    ap.add_argument("--cross-only", dest="bCrossOnly", action="store_true", 
                    help="Draw only cross-hemisphere connections (for hemi == bh or xh only)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)

    # === Parse input arguments === #
    args = ap.parse_args()

    componentFile = args.componentFile
    outputImgFN = args.outputImgFN
    hemi = args.hemi
    t_opacity = args.t_opacity
    bNoText = args.bNoText

    if args.bCrossOnly and not (hemi == "bh" or hemi == "xh"):
        error_log("--cross-only used with hemi other than bh or xh")
    
    if not (hemi == "lh" or hemi == "rh" or hemi == "bh" or hemi == "xh"):
        raise Exception, "Unexpected hemi: %s" % hemi

    if args.altStructVol != None:
        STRUCT_VOL = args.altStructVol

    if args.altCoordFile != None:
        COORD_FILE = args.altCoordFile

    # print(hemi) # DEBUG

    # === Load the structural image === #
    check_file(STRUCT_VOL)
    sImg = nb.load(STRUCT_VOL)
    sImgDat = sImg.get_data()

    # DEBUG
    # sImgDat = sImgDat[0 : 4 : 256, 0 : 4: 256, 0 : 4: 256]

    # === Downsample === #
    D = 2
    N = len(sImgDat) / D
    sImgDat1 = np.zeros([N, N, N])

    for i0 in range(0, N):
        for i1 in range(0, N):
            for i2 in range(0, N):
                sImgDat1[i0, i1, i2] = sImgDat[i0 * D, i1 * D, i2 * D]
            
    sImgDat = sImgDat1;
    
    #sInds = np.mgrid[0 : sImg.shape[0], 0 : sImg.shape[1], 0 : sImg.shape[2]]
    sInds = np.mgrid[0 : len(sImgDat), 
                     0 : len(sImgDat[0]), 
                     0 : len(sImgDat[0][0])]

    imgMin = np.min(sImgDat)
    imgMax = np.max(sImgDat)
    print("Image min = %f" % imgMin)
    print("Image min = %f" % imgMax)

    # mlab.figure(size=(600, 450))


    # === Read coordinates === #
    check_file(COORD_FILE)

    cf = open(COORD_FILE, "rt")
    ct = cf.read().split('\n')
    ct = remove_empty_strings(ct)
    cf.close()

    roi_names = []
    roi_nums = []
    roi_coords = []
    for (i0, tline) in enumerate(ct):
        t_items = tline.split(' ')
        if len(t_items) != 5:
            raise Exception, "Unrecognized formant in a line of %s: %s" \
                             % (COORD_FILE, tline)
        roi_names.append(t_items[0])
        roi_nums.append(t_items[1])
        t_coord = [float(t_items[2]), float(t_items[3]), float(t_items[4])]
        roi_coords.append(t_coord)

    assert(len(roi_names) == len(roi_nums))
    assert(len(roi_names) == len(roi_coords))

    # === Read component file === #
    check_file(componentFile)
    compf = open(componentFile, "rt")
    compt = compf.read()
    compf.close()
    compt = remove_empty_strings(compt.split('\n'))

    edges = []
    compNums = [] # Component number
    sigVals = []
    for (i0, tline) in enumerate(compt):
        assert(tline.count(" - ") == 1)
        assert(tline.count(": sig=") == 1)

        linkName = tline.split(": sig=")[0]

        roi1Name = linkName.split(" - ")[0]
        roi2Name = linkName.split(" - ")[1]

        if hemi == "xh" or hemi == "bh":
            hemi1 = roi1Name.split("_")[0]
            hemi2 = roi2Name.split("_")[0]
            if args.bCrossOnly and (hemi1 == hemi2):
                print("Skipping same-hemisphere projection: %s --> %s" % 
                      (roi1Name, roi2Name))
                continue

        if tline.count(", ") == 1:
            compNums.append(int(tline.split(", ")[1]))
            tline = tline.split(", ")[0]
        else:
            compNums.append(1)

        sigVal = np.abs(float(tline.split(": sig=")[1]))
                
        if args.bTranslateROIName:
            roi1Name = translate_roi_name(roi1Name)
            roi2Name = translate_roi_name(roi2Name)
        
        roi1Num = roi_names.index(roi1Name)
        roi2Num = roi_names.index(roi2Name)
        
        edges.append([roi1Num, roi2Num])
        sigVals.append(sigVal)

    assert(len(edges) == len(sigVals))
    assert(len(edges) == len(compNums))

    # === Call mayavi for 3d drawing === #
    # edges = edges[:4] # DEBUG 

    labPlotted = [0] * len(roi_names)
    
    # mlab.plot3d([0, 100], [0, 100], [0, 100], tube_radius=2.5)   
    # === Plot the "axes" === #
    """
    mlab.plot3d([-128, 128], [-128, -128], [-128, -128], \
                color=(1, 0, 0), tube_radius=2)   
    mlab.plot3d([-128, 128], [128, 128], [-128, -128], \
                color=(1, 1, 0), tube_radius=2)
    mlab.plot3d([-128, 128], [-128, -128], [128, 128], \
                color=(0, 1, 0), tube_radius=2)
    mlab.plot3d([-128, 128], [128, 128], [128, 128], \
                color=(0, 1, 0), tube_radius=2)

    mlab.plot3d([-128, -128], [-128, 128], [-128, -128], \
                color=(0, 0, 1), tube_radius=2)
    mlab.plot3d([128, 128], [-128, 128], [-128, -128], \
                color=(0, 1, 1), tube_radius=2)
    """

    mlab.figure(size=(800, 600), bgcolor=(1.0, 1.0, 1.0))

    # Coordinates: (in increasing coordinate value)
    #    Dimension 1: left to right (In image: R to L: Need - )
    #    Dimension 2: superior to inferior (In image: S to I: Okay )
    #    Dimension 3: anterior to posterior (In image: P to A: Need - )
    #  Green: anterior
    #  Red: superior
    #  Blue: right

    # mlab.show()
    
    #"""
    for (i0, tlink) in enumerate(edges):
        t_x = np.array([roi_coords[tlink[0]][0], \
                        roi_coords[tlink[1]][0]]) / D - N / 2
        t_x = -t_x
        t_y = np.array([roi_coords[tlink[0]][1], \
                        roi_coords[tlink[1]][1]]) / D- N / 2
        t_z = np.array([roi_coords[tlink[0]][2], \
                        roi_coords[tlink[1]][2]]) / D - N / 2

        mlab.plot3d(t_x, t_y, t_z, tube_radius=sigVals[i0] * 0.3 / D, \
                    color=COMPONENT_CLRS[compNums[i0] - 1], \
                    opacity=t_opacity)
        # mlab.plot3d(t_x, t_y, t_z, tube_radius=1)

        for j in range(2):
            if labPlotted[tlink[j]] == 0:
                mlab.points3d(t_x[j], t_y[j], t_z[j], \
                              scale_factor = 2.0 / D, \
                              color=COMPONENT_CLRS[compNums[i0] - 1], \
                              opacity=t_opacity)

                if hemi == "lh" or hemi == "rh":
                    rName = roi_names[tlink[j]].replace("lh_", "")\
                                               .replace("rh_", "")
                else:
                    rName = roi_names[tlink[j]].replace("lh_", "L ")\
                                               .replace("rh_", "R ")
                

                if not bNoText:
                    mlab.text3d(t_x[j], t_y[j], t_z[j], \
                                rName, \
                                color=(0, 0, 0), scale=2.0 / D)
                labPlotted[tlink[j]] = 1
    #"""

#        points3d(t_x, t_y, t_z)

#    mlab.roll(-90)
#    mlab.pitch(90)
#    mlab.move([600, 600, 400])

    #"""
    src = mlab.pipeline.scalar_field((sInds[0] - N / 2), \
                                     (sInds[1] - N / 2), \
                                     (sInds[2] - N / 2), \
                                     sImgDat)
#    mlab.pipeline.volume(src, vmin=10, vmax=1600, 
#                         color=(1.0, 1.0, 1.0))
    mlab.pipeline.volume(src, vmin=10, vmax=args.vmax, 
                         color=(1.0, 1.0, 1.0))
    #"""
              
    
    if hemi == "lh":
        mlab.view(azimuth=210, elevation=100, roll=180, \
                  focalpoint=[0, 0, 0], distance=220 / D)
    elif hemi == "rh":
        mlab.view(azimuth=145, elevation=-80, roll=180, \
                  focalpoint=[0, 0, 0], distance=220 / D)
    elif hemi == "xh" or hemi == "bh":
        mlab.view(azimuth=90, elevation=10, roll=0, \
                  focalpoint=[0, 0, 0], distance=220 / D)
    
    #cam,foc = mlab.move()
    #print(cam)
    #print(foc)
    
#    view=mlab.view()
#    print(view)
    
    print("Rendering done.")

    # outputImgFN = "netw_component.png"
    os.system("rm -f %s" % outputImgFN)

    #"""
    print("Saving image file ...")
    mlab.savefig(outputImgFN)
    check_file(outputImgFN)
    print("Image file saved at: %s" % outputImgFN)
    #"""

    mlab.show()

    
