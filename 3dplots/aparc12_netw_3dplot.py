#!/usr/bin/python

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

COMPONENT_CLRS = [(1.0, 1.0, 1.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.25)]

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

    if len(sys.argv) == 1:
        ap.print_help()

    # === Parse input arguments === #
    args = ap.parse_args()
    componentFile = args.componentFile
    outputImgFN = args.outputImgFN
    hemi = args.hemi
    t_opacity = args.t_opacity
    bNoText = args.bNoText
    
    if not (hemi == "lh" or hemi == "rh"):
        raise Exception, "Unexpected hemi: %s" % hemi

    # print(hemi) # DEBUG

    # === Load the structural image === #
    check_file(STRUCT_VOL)
    sImg = nb.load(STRUCT_VOL)
    sImgDat = sImg.get_data()
    
    sInds = np.mgrid[0 : sImg.shape[0], 0 : sImg.shape[1], 0 : sImg.shape[2]]

    imgMin = np.min(sImgDat)
    imgMax = np.max(sImgDat)
    print(imgMin)
    print(imgMax)

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
        if tline.count(", ") == 1:
            compNums.append(int(tline.split(", ")[1]))
            tline = tline.split(", ")[0]
        else:
            compNums.append(1)

        assert(tline.count(" - ") == 1)
        assert(tline.count(": sig=") == 1)
        
        linkName = tline.split(": sig=")[0]
        sigVal = np.abs(float(tline.split(": sig=")[1]))

        roi1Name = linkName.split(" - ")[0]
        roi2Name = linkName.split(" - ")[1]
        
        roi1Num = roi_names.index(roi1Name)
        roi2Num = roi_names.index(roi2Name)
        
        edges.append([roi1Num, roi2Num])
        sigVals.append(sigVal)

    assert(len(edges) == len(sigVals))
    assert(len(edges) == len(compNums))

    # === Call mayavi for 3d drawing === #
    # edges = edges[:4] # DEBUG 

    labPlotted = [0] * len(roi_names)

#    sImgDat = sImgDat[::-1, :, :]

    
    
        
    
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
    # sys.exit(0)

    for (i0, tlink) in enumerate(edges):
        t_x = np.array([roi_coords[tlink[0]][0], \
                        roi_coords[tlink[1]][0]]) - 128
        t_x = -t_x
        t_y = np.array([roi_coords[tlink[0]][1], \
                        roi_coords[tlink[1]][1]]) - 128
        t_z = np.array([roi_coords[tlink[0]][2], \
                        roi_coords[tlink[1]][2]]) - 128

        mlab.plot3d(t_x, t_y, t_z, tube_radius=sigVals[i0] * 0.3, \
                    color=COMPONENT_CLRS[compNums[i0] - 1], \
                    opacity=t_opacity)
        # mlab.plot3d(t_x, t_y, t_z, tube_radius=1)

        for j in range(2):
            if labPlotted[tlink[j]] == 0:
                mlab.points3d(t_x[j], t_y[j], t_z[j], \
                              scale_factor=2.0, \
                              color=COMPONENT_CLRS[compNums[i0] - 1], \
                              opacity=t_opacity)
                if not bNoText:
                    mlab.text3d(t_x[j], t_y[j], t_z[j], \
                                roi_names[tlink[j]].replace("lh_", "").replace("rh_", ""), \
                                color=(0, 0, 0), scale=2.0)
                labPlotted[tlink[j]] = 1

#        points3d(t_x, t_y, t_z)

#    mlab.roll(-90)
#    mlab.pitch(90)
#    mlab.move([600, 600, 400])

    #"""
    src = mlab.pipeline.scalar_field((sInds[0] - 128), \
                                     (sInds[1] - 128), \
                                     (sInds[2] - 128), \
                                     sImgDat)
    mlab.pipeline.volume(src, vmin=10, vmax=400) # vmax=500 is good?
    #"""
              
    if hemi == "lh":
        mlab.view(azimuth=200, elevation=100, roll=180, \
                  focalpoint=[0, 0, 0], distance=250)
    else:
        mlab.view(azimuth=160, elevation=-80, roll=180, \
                  focalpoint=[0, 0, 0], distance=250)
    



    cam,foc = mlab.move()
    print(cam)
    print(foc)
    
    view=mlab.view()
    print(view)
    
    # outputImgFN = "netw_component.png"
    os.system("rm -f %s" % outputImgFN)
    print("Saving image file ...")
    mlab.savefig(outputImgFN)
    check_file(outputImgFN)
    print("Image file saved at: %s" % outputImgFN)
    mlab.show()

    
 
  

    sys.exit(0)
    
