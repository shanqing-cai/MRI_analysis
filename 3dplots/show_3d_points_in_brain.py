#!/usr/bin/python

import os
import sys
import numpy as np
import nibabel as nb
import argparse
import mayavi.mlab as mlab

from scai_utils import *

COORD_FILE = "/home/cais/STUT/FSDATA/fsaverage2/mri/aparc12_roi_coords.txt"
STRUCT_VOL = "/home/cais/STUT/FSDATA/MNI152_1mm/mri/brain.nii.gz"
# STRUCT_VOL = "/home/cais/STUT/FSDATA/MNI152_1mm/mri/brain.nii.gz"

# --- To convert into input --- #
# componentFile = "corrSSI4_sigComponents_1.txt"
# edges = [[0, 1], [1, 2], [2, 3], [3, 0], [1, 3]]

mniCoord = [-41, 8, 16]
ptClr = (0, 1, 0)
#coord = [168, 131, 153]

AZIMUTH = {"lh": 210, "rh": 160}
ELEVATION = {"lh": 100, "rh": -80}

DEFAULT_ANGLE = -1000.0

DEFAULT_RADIUS = 3.0

# MNI-to-raster coordinate transform
def mni2raster(mniCoord):
    assert(len(mniCoord) == 3)

    coord = np.zeros(3)
    coord[0] = 127.0 - mniCoord[0]
    coord[1] = 147.0 - mniCoord[2]
    coord[2] = 145.0 + mniCoord[1]

    return coord

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Plot 3D points (as spheres, of course) inside a translucent template brain (STRUCT_VOL = %s)" % STRUCT_VOL)
    ap.add_argument("ptsListFN", help="Input points list file name")
    ap.add_argument("--azimuth", dest="arg_azimuth", \
                    type=float, default=DEFAULT_ANGLE, help="Azimuth angle")
    ap.add_argument("--elevation", dest="arg_elevation", \
                    type=float, default=DEFAULT_ANGLE, help="Azimuth angle")
    ap.add_argument("--radius", dest="radius", type=float, \
                    default=DEFAULT_RADIUS, \
                    help="Radius of the spheres (default: %.1f)" \
                         % DEFAULT_RADIUS)
    ap.add_argument("--noBrain", dest="bNoBrain", action="store_true", \
                    help="Do not show the brain (usually used for speed)")
    ap.add_argument("--out", type=str, dest="outputImgFN", default="")

#    ap.add_argument("componentFile", type=str)
#    ap.add_argument("outputImgFN", type=str)
    ap.add_argument("--hemi", type=str, default="lh")


    # === Parse input arguments === #
    args = ap.parse_args()
    ptsListFN = args.ptsListFN
    outputImgFN = args.outputImgFN
    hemi = args.hemi

    arg_azimuth = args.arg_azimuth
    arg_elevation = args.arg_elevation

    radius = args.radius
    
    bNoBrain = args.bNoBrain

    # === Load input points (MNI coordinates and colors) === #
    check_file(ptsListFN)

    ptsMNIs = []
    ptsClrs = []
    ptsListT = remove_empty_strings(read_text_file(ptsListFN))
    
    for (i0, tline) in enumerate(ptsListT):
        if tline.startswith("#"):
            continue

        t_items = remove_empty_strings(tline.replace('\t', ' ').split(' '))
        assert(len(t_items) >= 6)

        ptsMNIs.append([float(t_items[0]), \
                        float(t_items[1]), \
                        float(t_items[2])])
        ptsClrs.append(tuple([float(t_items[3]), \
                              float(t_items[4]), \
                              float(t_items[5])]))

    # === Load the structural image === #
    check_file(STRUCT_VOL)
    sImg = nb.load(STRUCT_VOL)
    sImgDat = sImg.get_data()
    
    sInds = np.mgrid[0 : sImg.shape[0], 0 : sImg.shape[1], 0 : sImg.shape[2]]

    imgMin = np.min(sImgDat)
    imgMax = np.max(sImgDat)
    print(imgMin)
    print(imgMax)

    # === Read coordinates === #
    
    mlab.figure(size=(800, 600), bgcolor=(1.0, 1.0, 1.0))
#    mlab.figure(size=(800, 600), bgcolor=(0.0, 0.0, 0.0))

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

    # Coordinates: (in increasing coordinate value)
    #    Dimension 1: left to right (In image: R to L: Need - )
    #    Dimension 2: superior to inferior (In image: S to I: Okay )
    #    Dimension 3: anterior to posterior (In image: P to A: Need - )
    #  Green: anterior
    #  Red: superior
    #  Blue: right

    for (i0, ptMNI) in enumerate(ptsMNIs):
        coord = mni2raster(ptMNI)
        
        t_x = -(coord[0] - 128)
        t_y = coord[1] - 128 
        t_z = coord[2] - 128

        print("Plotting point: [%.1f, %.1f, %.1f]" % (t_x, t_y, t_z))

        mlab.points3d(t_x, t_y, t_z, scale_factor=radius, color=ptsClrs[i0])
#    mlab.points3d(0, 0, 0, scale_factor=10.0, colormap="Reds")

    #"""
    if not bNoBrain:
        src = mlab.pipeline.scalar_field((sInds[0] - 128), \
                                         (sInds[1] - 128), \
                                         (sInds[2] - 128), \
                                         sImgDat)
        mlab.pipeline.volume(src, vmin=10, vmax=2500, 
                             color=(1.00, 1.00, 1.00)) # vmax=500 is good?
    #"""

    if arg_azimuth != DEFAULT_ANGLE:
        t_azimuth = arg_azimuth
    else:
        t_azimuth = AZIMUTH[hemi]

    if arg_elevation != DEFAULT_ANGLE:
        t_elevation = arg_elevation
    else:
        t_elevation = ELEVATION[hemi]

    print("azimuth = %f; elevation = %f" % (t_azimuth, t_elevation))
    
    mlab.view(azimuth=t_azimuth, elevation=t_elevation, roll=180, \
              focalpoint=[0, 0, 0], distance=240)


    cam,foc = mlab.move()
    print(cam)
    print(foc)
    
    view=mlab.view()
    print(view)

    # === Save to output image === #
    if len(outputImgFN) > 0:
        os.system("rm -f %s" % outputImgFN)
        print("Saving image file ...")
        mlab.savefig(outputImgFN)
        check_file(outputImgFN)
        print("Image file saved at: %s" % outputImgFN)
    
    mlab.show()
    
