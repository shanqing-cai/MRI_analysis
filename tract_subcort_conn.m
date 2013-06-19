function tract_subcort_conn(scSeed, mainMaskRatio)
%% Constants
tractSubcortConnDir = '/users/cais/STUT/analysis/tract_subcort_conn';

grps = {'PWS', 'PFS'};

%% Config
P_THRESH_UNC = 0.05;

%%
check_dir(tractSubcortConnDir);

seedDir = fullfile(tractSubcortConnDir, scSeed);
check_dir(seedDir);

% -- Load the list of cortical ROIs -- %
rois = aparc12_cortical_rois();

%% -- Load data -- %%
sccConn = struct;
sccConn.PWS = [];
sccConn.PFS = [];

for i1 = 1 : numel(rois)
    roi = rois{i1};
    roiDir = fullfile(seedDir, roi);
    
    if isdir(roiDir)
        stfn = fullfile(roiDir, sprintf('stats_%.2f.txt', mainMaskRatio));
        check_file(stfn);
        stt = textread(stfn, '%s', 'delimiter', '\n');
        
        for k1 = 1 : numel(grps)
            grp = grps{k1};
            
            assert(isequal(stt{2 + (k1 - 1) * 2}(1 : 3), grp));
            nums = splitstring(stt{3 + (k1 - 1) * 2});
            if size(sccConn.(grp), 2) == 0
                sccConn.(grp) = nan(length(nums), length(rois));            
            end
        
            for i2 = 1 : numel(nums)
                sccConn.(grp)(i2, i1) = str2double(nums{i2});
            end                   
        end
    else
        fprintf(2, 'WARNING: Cannot find directory for ROI %s\n', roi);
    end
end

%% Statistical comparisons
rs_zs = nan(1, numel(rois));
rs_ps = nan(1, numel(rois));
rs_rs = nan(1, numel(rois));
rs_sigs = nan(1, numel(rois));

medDiffs = nan(1, numel(rois));

hiliteROIs = {};
hiliteROIClrs = {};
drawBndClrs = {};

for i1 = 1 : numel(rois)
    roi = rois{i1};
    
    [rs_ps(i1), ~, stats] = ranksum(sccConn.PWS(:, i1), sccConn.PFS(:, i1));
    rs_zs(i1) = stats.zval;
    rs_rs(i1) = stats.ranksum;
    
    medDiffs(i1) = median(sccConn.PWS(:, i1)) - median(sccConn.PFS(:, i1));
    rs_sigs(i1) = sign(medDiffs(i1)) * -log10(rs_ps(i1));
    
    if rs_ps(i1) < P_THRESH_UNC
        fprintf(1, '%s: ranksum z=%.3f; p=%.3f\n\t(Median: PWS=%f; PFS=%f)\n\t(IQR; PWS=%f; PFS=%f)\n', ...
                roi, rs_zs(i1), rs_ps(i1), ...
                median(sccConn.PWS(:, i1)), median(sccConn.PFS(:, i1)), ...
                iqr(sccConn.PWS(:, i1)), iqr(sccConn.PFS(:, i1)));
            
        hiliteROIs{end + 1} = roi;
        hiliteROIClrs{end + 1} = [0, 1, 0];
        drawBndClrs{end + 1} = [0, 1, 0];
    end
    
%     t_fillROIs{end + 1} = roi;
%     t_fillClrs{end + 1} = 
end

%% Visualization
cm0 = create_green_red_colormap();
cm = nan(size(cm0, 1) + 1, size(cm0, 2));
for i1 = 1 : size(cm, 2)
    cm(:, i1) = interp1(1 : size(cm0, 1), cm0(:, i1), ...
                        linspace(1, size(cm0, 1), size(cm, 1)));
end

t_fillROIs = rois;
t_fillClrs = {};



max_abs_sig = nanmax(abs(rs_sigs));
for i1 = 1 : numel(rois)
    if ~isnan(rs_sigs(i1))
        tpos = (rs_sigs(i1) + max_abs_sig) / 2 / max_abs_sig;
        tpos = round(1 + tpos * (size(cm, 1) - 1));       
        t_fillClrs{i1} = cm(tpos, :);
    else
        t_fillClrs{i1} = [1, 1, 1];
    end
end

hemi = scSeed(1 : 2);
roiFigs.(hemi) = ...
        sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc_noText.tif', hemi);
roiFigsDiv.(hemi) = ...
        sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc_noText_div.tif', hemi);
roiFigsText.(hemi) = ...
        sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc.tif', hemi);
    
check_file(roiFigs.(hemi));
check_file(roiFigsDiv.(hemi));
check_file(roiFigsText.(hemi));

imBase = imread(roiFigs.(hemi));
imDiv = imread(roiFigsDiv.(hemi));
imText = imread(roiFigsText.(hemi));

imFASigMap.(hemi) = ...
        fill_aparc12_roi_fig(imBase, imDiv, imText, hemi, ...
                             t_fillROIs, 'fillClrs', t_fillClrs, ...
                             'clrNameROIs', hiliteROIs, hiliteROIClrs, ...
                             'drawBndROIs', hiliteROIs, drawBndClrs);

figure;
imshow(imFASigMap.(hemi));

return