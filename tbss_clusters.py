#!/usr/bin/python

import os
import sys
import glob
import argparse
import tempfile
import numpy as np
from scipy import stats
from subprocess import Popen, PIPE
import xml.etree.ElementTree as ET

from scai_utils import *
from get_qdec_info import get_qdec_info

atlas_label_fn = \
    "/usr/share/fsl/5.0/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz"
atlas_tract_fn = \
    "/usr/share/fsl/5.0/data/atlases/JHU/JHU-ICBM-tracts-prob-1mm.nii.gz"
atlas_label_xml = \
    "/usr/share/fsl/5.0/data/atlases/JHU-labels.xml"
atlas_tract_xml = \
    "/usr/share/fsl/5.0/data/atlases/JHU-tracts.xml"

aparc12_full_ctab = "/users/cais/STUT/scripts/slaparc_550.ctab"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Get cluster summary from TBSS t-statistic files (*_tstat?.nii.gz)")
    ap.add_argument("tstatfn", type=str, \
                        help="tstat image file (.nii.gz format)")
    ap.add_argument("voxp", type=float, \
        help="Voxel-wise p-value threshold, two-tailed (e.g., 0.001)")
    ap.add_argument("voxcnt", type=int, \
        help="Voxel counter threshold (Unit: ) (e.g., 10)")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    tstatfn = args.tstatfn
    voxp = args.voxp
    voxcnt = args.voxcnt

    # Input sanity check
    if not tstatfn.startswith("/"):
        tstatfn = os.path.abspath(tstatfn)

    if voxp <= 0 or voxp >= 1:
        raise Exception, "Invalid value of voxp: %f" % voxp
    
    if voxcnt <= 0:
        raise Exception, "Invalid value of voxcnt: %d" % voxcnt

    check_file(tstatfn)

    # === Read xml files for labels == # 
    check_file(atlas_label_xml)
    tree = ET.parse(atlas_label_xml)
    labs = tree.getroot()
    a = labs[1]
    b = a.getchildren()

    atl_labs = {"ind": [], "name": []}
    for tb in b:
        atl_labs["ind"].append(int(tb.attrib["index"]))
        atl_labs["name"].append(tb.text)

    # === Read xml files for tracts === # 
    check_file(atlas_tract_xml)
    tree = ET.parse(atlas_tract_xml)
    labs = tree.getroot()
    a = labs[1]
    b = a.getchildren()

    atl_tracts = {"ind": [], "name": []}
    for tb in b:
        atl_tracts["ind"].append(int(tb.attrib["index"]))
        atl_tracts["name"].append(tb.text)

    # === Read full aparc12 (SLaparc) color table === %
    check_file(aparc12_full_ctab)
    (roi_nums, roi_names) = read_ctab(aparc12_full_ctab)
    
    #sys.exit(0)
    
    # Search for the all_FA.nii.gz file, for determining the number of subjects
    (fpath, fn) = os.path.split(tstatfn)
    allFA = os.path.join(fpath, "all_FA.nii.gz")
    check_file(allFA)
    (sout, serr) = Popen(["mri_info", allFA], \
                         stdout=PIPE, stderr=PIPE).communicate()
    if len(serr) > 0:
        raise Exception, "ERROR occurred during mri_info %s" % allFA
    sout = sout.split('\n')
    N = int(sout[2].split(' ')[-1])

    df = N - 2
    print("INFO: N = %d; df = %d" % (N, df))

    # Binarize
    tthr = -stats.t.ppf(voxp / 2.0, df)
    print("INFO: t-value thr = %f" %  tthr)

    bin_out = os.path.join(fpath, fn.replace(".nii.gz", \
                                             "_pthr%f.nii.gz" % voxp))
    bin_cmd = "mri_binarize --i %s --min %f --o %s" % \
              (tstatfn, tthr, bin_out)
    os.system("rm -f %s" % bin_out)
    saydo(bin_cmd)
    check_file(bin_out)

    # Get masked t-value file
    masked_fn = os.path.join(fpath, fn.replace(".nii.gz", \
                                               "_pthr%f.masked.nii.gz" % voxp))
    os.system("rm -f %s" % masked_fn)
    mul_cmd = "fslmaths %s -mul %s %s" % (tstatfn, bin_out, masked_fn)
    saydo(mul_cmd)
    check_file(masked_fn)
    
    # Run mri_volcluster
    sum_fn = os.path.join(fpath, fn.replace(".nii.gz", \
                                 "_pthr%f_cnt%d.sum" % (voxp, voxcnt)))
    volclust_out = os.path.join(fpath, fn.replace(".nii.gz", \
                                 "_pthr%f_cnt%d.vc.nii.gz" \
                                 % (voxp, voxcnt)))
    mvc_cmd = "mri_volcluster --in %s --thmin 0.5 --minsizevox %d --sum %s --out %s" % \
              (masked_fn, voxcnt, sum_fn, volclust_out)
    print("sum_fn = %s" % sum_fn)
#    os.system("rm -f %s" % sum_fn)
    os.system("rm -f %s" % volclust_out)
    saydo(mvc_cmd)
    check_file(sum_fn)
    check_file(volclust_out)

    # Get vc-masked t-value file
    vcmasked_fn = os.path.join(fpath, \
                               fn.replace(".nii.gz", \
                               "_pthr%f.vcmasked.nii.gz" % voxp))
    os.system("rm -f %s" % vcmasked_fn)
    mul_cmd = "fslmaths %s -mul %s %s" % (masked_fn, volclust_out, vcmasked_fn)
    saydo(mul_cmd)
    check_file(vcmasked_fn)
    
    # === Load the summary file === #
    sum_f = open(sum_fn, "r")
    sumt = sum_f.read().split('\n')
    sum_f.close()
    sumt = remove_empty_strings(sumt)
    
    nClust = 0
    clustSizes = []
    clustSizesVox = []
    clustX = []
    clustY = []
    clustZ = []
    clust_mniX = []
    clust_mniY = []
    clust_mniZ = []
    maxT = []
    maxCohenD = []

    for (i0, tline) in enumerate(sumt):
        if tline[0] == "#":
            continue
        
        t_items = tline.replace('\t', ' ').split(' ')
        t_items = remove_empty_strings(t_items)

        if len(t_items) != 7:
            raise Exception, "Unrecognized format in line: %s" % tline

        nClust = nClust + 1
        clustSizesVox.append(int(t_items[1]))
        clustSizes.append(float(t_items[2]))

        clustX.append(float(t_items[3]))
        clustY.append(float(t_items[4]))
        clustZ.append(float(t_items[5]))

        clust_mniX.append(90.0 - clustX[-1] * 1.0)
        clust_mniY.append(-126.0 + clustY[-1] * 1.0)
        clust_mniZ.append(-72.0 + clustZ[-1] * 1.0)

        maxT.append(float(t_items[6]))

    # === Determine the labels and tracts of the clusters === #
    clustLabels = [""] * nClust
    clustTracts = [""] * nClust
    clustAparc12Lab = [""] * nClust

    (tbssDir, foo) = os.path.split(tstatfn)
    (tbssDir, foo) = os.path.split(tbssDir)
    mergedAparc12Lab = os.path.join(tbssDir, "aparc12", "merged.nii.gz")
    if not os.path.isfile(mergedAparc12Lab):
        saydo("gen_tbss_aparc12_prob_map.py %s" % tbssDir)

    # === Locate the all_FA_skeletonised images (for calculating z-scores) === #
    aFASkel = os.path.join(tbssDir, "stats", "all_FA_skeletonised.nii.gz")
    check_file(aFASkel)

    # === Find out the subject IDs and their groups === #
    origDir = os.path.join(tbssDir, "origdata")
    check_dir(origDir)

    ds = glob.glob(os.path.join(origDir, "S??.nii.gz"))
    ds.sort()

    sIDs = []
    idxPWS = []
    idxPFS = []

    for (i0, d) in enumerate(ds):
        [tpath, tfn] = os.path.split(d)
        sID = tfn.replace(".nii.gz", "")
        sIDs.append(sID)
        
        if get_qdec_info(sID, "diagnosis") == "PWS":
            idxPWS.append(i0)
        elif get_qdec_info(sID, "diagnosis") == "PFS":
            idxPFS.append(i0)
        else:
            raise Exception, "Unrecognized diagnosis for subject %s: %s" % \
                             (sID, get_qdec_info(sID, "diagnosis"))
    
    # === Process the clusters === #
    for i0 in range(nClust):
        # == Determine label == #
        roi_fn = tempfile.mktemp() + ".nii.gz"
        roi_cmd = "fslroi %s %s %d 1 %d 1 %d 1" % \
                  (atlas_label_fn, roi_fn, clustX[i0], clustY[i0], clustZ[i0])
        saydo(roi_cmd)
        check_file(roi_fn)

        (sout, serr) = Popen(["fslstats", roi_fn, "-M"], \
                             stdout=PIPE, stderr=PIPE).communicate()
        if len(serr) > 0:
            raise Exception, "ERROR occurred during fslstats on file %s" % \
                             roi_fn
        
        labn = int(np.round(float(sout.split(' ')[0])))

        clustLabels[i0] = atl_labs['name'][atl_labs['ind'].index(labn)]

        os.system("rm -f %s" % roi_fn)
        
        # == Determine tract == #
        roi_fn = tempfile.mktemp() + ".nii.gz"
        roi_cmd = "fslroi %s %s %d 1 %d 1 %d 1" % \
                  (atlas_tract_fn, roi_fn, clustX[i0], clustY[i0], clustZ[i0])
        saydo(roi_cmd)
        check_file(roi_fn)

        (sout, serr) = Popen(["fslstats", "-t", roi_fn, "-M"], \
                             stdout=PIPE, stderr=PIPE).communicate()
        if len(serr) > 0:
            raise Exception, "ERROR occurred during fslstats on file %s" % \
                             roi_fn

        items = sout.replace('\n', ' ').split(' ')
        items = remove_empty_strings(items)
        vals = []
        for item in items:
            vals.append(float(item))

        if len(vals) != len(atl_tracts['ind']):
            raise Exception, "Unexpected number of frames in file: %s" % \
                             atlas_tract_fn
        
        if vals.count(0.0) == len(vals):
            clustTracts[i0] = "Undetermined"
        else:
            t_max = np.max(vals)
            t_idx = vals.index(t_max)
            
            clustTracts[i0] = atl_tracts['name'][atl_tracts['ind'].index(t_idx)]

        os.system("rm -f %s" % roi_fn)

        # == Determine dominant aparc12 label == #
        roi_fn = tempfile.mktemp() + ".nii.gz"
        roi_cmd = "fslroi %s %s %d 1 %d 1 %d 1" % \
                  (mergedAparc12Lab, roi_fn, clustX[i0], clustY[i0], clustZ[i0])
        saydo(roi_cmd)
        check_file(roi_fn)

        (sout, serr) = Popen(["fslstats", "-t", roi_fn, "-M"], \
                             stdout=PIPE, stderr=PIPE).communicate()

        items = sout.replace('\n', ' ').split(' ')
        items = remove_empty_strings(items)
        vals = []
        for item in items:
            vals.append(float(item))

        counts = np.bincount(vals)
        idxmax = np.argmax(counts)

        if roi_nums.count(idxmax) == 1:
            clustAparc12Lab[i0] = roi_names[roi_nums.index(idxmax)]
        else:
            clustAparc12Lab[i0] = np.nan

        os.system("rm -f %s" % roi_fn)

        # == Determine the z-value == #
        roi_fn = tempfile.mktemp() + ".nii.gz"
        roi_cmd = "fslroi %s %s %d 1 %d 1 %d 1" % \
                  (aFASkel, roi_fn, clustX[i0], clustY[i0], clustZ[i0])
        saydo(roi_cmd)
        check_file(roi_fn)

        (sout, serr) = Popen(["fslstats", "-t", roi_fn, "-M"], \
                             stdout=PIPE, stderr=PIPE).communicate()

        items = sout.replace('\n', ' ').split(' ')
        items = remove_empty_strings(items)
        vals = []
        for item in items:
            vals.append(float(item))

        vals = np.array(vals)    
        vals_PWS = vals[idxPWS]
        vals_PFS = vals[idxPFS]
        mean_PWS = np.mean(vals_PWS)
        mean_PFS = np.mean(vals_PFS)
        std_PWS = np.std(vals_PWS)
        std_PFS = np.std(vals_PFS)
        std_2g = np.sqrt(((len(vals_PWS) - 1) * std_PWS * std_PWS + \
                          (len(vals_PFS) - 1) * std_PFS * std_PFS) / \
                             (len(vals_PWS) + len(vals_PFS) - 2))
        maxCohenD.append((mean_PWS - mean_PFS) / std_2g)

        os.system("rm -f %s" % roi_fn)
    
    
    # --- Print viewing command --- #
    mean_FA = os.path.join(fpath, "mean_FA.nii.gz")
    check_file(mean_FA)

    skel_mask = os.path.join(fpath, "mean_FA_skeleton_mask.nii.gz")
    check_file(skel_mask)

    check_file(atlas_label_fn)

    viewCmd = "freeview %s %s:colormap=nih %s:colormap=jet %s:colormap=nih:opacity=0.25" % \
              (mean_FA, skel_mask, vcmasked_fn, atlas_label_fn)
    
    print("------------------------------------------")
    print("\nTo view the results, do:\n\t%s" % viewCmd)

    print("\n")
    for i0 in range(nClust):
        print("Clust #%d:\n\tVolume coord = [%d, %d, %d]\n\tMNI coord = [%.1f, %.1f, %.1f]\n\tsize = %d voxels\n\tPeak t = %f\n\tPeak Cohen's d = %f\n\tlabel = %s\n\tMax tract = %s\n\tMax aparc12 label = %s" % \
              (i0 + 1, clustX[i0], clustY[i0], clustZ[i0], \
               clust_mniX[i0], clust_mniY[i0], clust_mniZ[i0], \
               clustSizes[i0], maxT[i0], maxCohenD[i0], \
               clustLabels[i0], clustTracts[i0], clustAparc12Lab[i0]))

    if nClust == 0:
        print("nClust = 0: Did not find any significant clusters")
