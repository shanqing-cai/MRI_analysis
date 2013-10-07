def get_parc_stats(fsDir, subjID, parcName, bVerbose=False):
    import sys, os
    import tempfile
    from scai_utils import check_dir, check_file, check_bin_path, \
                           info_log, cmd_stdout, saydo, read_text_file, \
                           remove_empty_strings

    #=== Constants ===#
    hemis = ["lh", "rh"]
    segStatsBin = "mri_segstats"

    #=== Implementation ===#
    check_bin_path(segStatsBin)
    
    sDir = os.path.join(fsDir, subjID)
    check_dir(sDir)

    lblDir = os.path.join(sDir, "label")
    check_dir(sDir)

    if bVerbose:
        info_log("Label directory = %s" % lblDir)

    morphInfo = {}

    rois = []
    area_mm2 = []
    
    for (i0, hemi) in enumerate(hemis):
        if bVerbose:
            info_log("Working on hemisphere: %s" % hemi)

        annot = os.path.join(lblDir, "%s.%s.annot" % (hemi, parcName))
        check_file(annot)
        
        tmpSum = tempfile.mktemp()
        cmd = "%s --annot %s %s %s --sum %s" % \
              (segStatsBin, subjID, hemi, parcName, tmpSum)
        saydo(cmd)
    
        check_file(tmpSum)    
        if bVerbose:
            print("Intermediate results saved at: %s" % tmpSum)

        t = read_text_file(tmpSum)

        for (i0, tline) in enumerate(t):
            tline = tline.strip()
            if len(tline) == 0:
                continue
            elif tline.startswith("#"):
                continue
            
            t_items = remove_empty_strings(tline.split())
            if len(t_items) != 5:
                continue

            t_roi = t_items[4]
            if t_roi.lower() == "none":
                continue

            rois.append("%s_%s" % (hemi, t_roi))
            area_mm2.append(float(t_items[3]))

        saydo("rm -rf %s" % tmpSum)
    morphInfo["rois"] = rois
    morphInfo["area_mm2"] = area_mm2

    return morphInfo
        
