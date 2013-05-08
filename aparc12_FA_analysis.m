function aparc12_FA_analysis(meas, xmm, dtiPrepMode)
sIDs.PWS = {'S01', 'S04', 'S06', 'S07', 'S08', 'S09', 'S10', 'S12', 'S15',  ...
            'S16', 'S17', 'S20', 'S21', 'S26', 'S28', 'S29', 'S33', 'S34', ...
            'S36', 'S37'};
            % Included S16 (But S16 should probably be kept in the group FA comparison) 
            % What's wrong with S21? Why was he excluded before, in
            % aparcSL_FA_analysis?
            % Left out S04: due to registration problems, and B0
            % distortions.
sIDs.PFS = {'S02', 'S03', 'S05', 'S11', 'S13', 'S14', 'S18', 'S19', 'S22', ...
            'S23', 'S25', 'S27', 'S30', 'S31', 'S32', 'S35', ...
            'S39'};
            % Left out S24 (gross structural brain abnormality)
 
% xmm = 3;
% ROI = 'lh_dIFo';

DATA_DIR = '/users/cais/STUT/DATA';

if isequal(dtiPrepMode, 'none')
    DAT_FN_WC = 'aparc12_%s_wm%dmm.mat';
elseif isequal(dtiPrepMode, 'dtiprep')
    DAT_FN_WC = 'aparc12_dtiprep_%s_wm%dmm.mat';
elseif isequal(dtiPrepMode, 'dtiprep2')
    DAT_FN_WC = 'aparc12_dtiprep2_%s_wm%dmm.mat';
else
    error('Unrecognized dtiPrepMode');
end

TTEST_P_THRESH_UC = 0.05;
LINCORR_P_THRESH_UC = 0.05;

%%
% ROIS_BEHAVCORR = {'lh_vIFo', 'rh_vIFo', ...
%                   'lh_dIFo', 'rh_dIFo', ...        
%                   'lh_vPMC', 'rh_vPMC', ...
%                   'lh_SMA', 'rh_SMA', ...
%                   'lh_preSMA', 'rh_preSMA', ...
%                   'lh_vMC', 'rh_vMC', ...
%                   'lh_vSC', 'rh_vSC', ...
%                   'lh_aCO', 'rh_aCO', ...
%                   'lh_pCO', 'rh_pCO', ...
%                   'lh_H',   'rh_H', ...
%                   'lh_pSTg', 'rh_pSTg', ...
%                   'lh_PT', 'rh_PT'};
ROIS_BEHAVCORR = {};

hemis = {'lh', 'rh'};

t_rois = get_aparc12_cortical_rois('speech');
for i1 = 1 : numel(hemis)
    hemi = hemis{i1};
    
    for i2 = 1 : numel(t_rois)
        ROIS_BEHAVCORR{end + 1} = sprintf('%s_%s', hemi, t_rois{i2});
    end
end

% ROIS_sel = {'lh.iFo', 'lh.vPMC', 'lh.vMC', 'lh.aCO', 'lh.vSSC', ...
%                 'lh.H', 'lh.PT'};
ROIS_sel = {'lh_vIFo', 'lh_vPMC', 'lh_vMC', 'lh_aCO', 'lh_vSC', ...
            'lh_H', 'lh_PT'};
% ROIS_sel = {'lh_aIFs', 'lh_vIFo', 'lh_vPMC', 'lh_vMC', 'lh_aCO', 'lh_vSC', ...
%             'lh_H', 'lh_PT'};
                            
%% Config: visualization
colors.PWS = 'r';
colors.PFS = 'k';

defaultFigClr = 'w';

figSaveDir = '/users/cais/STUT/figures';


for i0 = 1 : numel(hemis)
    hemi = hemis{i0};
    
    roiFigs.(hemi) = ...
        sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc_noText.tif', hemi);
    roiFigsDiv.(hemi) = ...
        sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc_noText_div.tif', hemi);
    roiFigsText.(hemi) = ...
        sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc.tif', hemi);
    
    check_file(roiFigs.(hemi));
    check_file(roiFigsDiv.(hemi));
    check_file(roiFigsText.(hemi));
end

hiliteROIClr = [0.5, 1, 0];
drawBndClr = [0, 1.0, 0];

%% Get aparc12 ROI names
%if nargin == 0
    roi_names_0 = get_aparc12_cortical_rois;
    roi_names = {};
    
    for hemi = {'lh', 'rh'};
        for i1 = 1 : numel(roi_names_0)
            roi_names{end + 1} = [hemi{1}, '_', roi_names_0{i1}];
        end
    end
%end

%% Get behavioral measures
grps = fields(sIDs);
SSI4 = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    auSTI.(grp) = get_qdec_measure(sIDs.(grp), 'auSTI');
    auSTI_n.(grp) = get_qdec_measure(sIDs.(grp), 'auSTI_n');
    auSTI_dsd.(grp) = get_qdec_measure(sIDs.(grp), 'auSTI_dsd');
    auSTI_dsd2.(grp) = get_qdec_measure(sIDs.(grp), 'auSTI_dsd2');
    hnSV.(grp) = get_qdec_measure(sIDs.(grp), 'hnSV');
    rnSV.(grp) = get_qdec_measure(sIDs.(grp), 'rnSV');
    nnSV.(grp) = get_qdec_measure(sIDs.(grp), 'nnSV');
           
    EHc.(grp) = get_qdec_measure(sIDs.(grp), 'EH_comp_300');
    tempResp.(grp) = get_qdec_measure(sIDs.(grp), 'tempResp');
    
    if isequal(grp, 'PWS')
        SSI4.(grp) = get_qdec_measure(sIDs.(grp), 'SSI');
    end
end

%%
nROIs = length(roi_names);
fa.PWS = nan(nROIs, numel(sIDs.PWS));
fa.PFS = nan(nROIs, numel(sIDs.PFS));

for h1 = 1 : nROIs
    t_roi = roi_names{h1};
    
    grps = fields(sIDs);
    for i1 = 1 : numel(grps)
        grp = grps{i1};

        for i2 = 1 : numel(sIDs.(grp))
            sID = sIDs.(grp){i2};
            matfn = fullfile(DATA_DIR, sID, sprintf(DAT_FN_WC, meas, xmm));
            if ~isfile(matfn)
                error('Cannot find mat file: %s', matfn)
            end
            load(matfn); % gives rois, meanfa

            iroi = strmatch(t_roi, rois, 'exact');
            if isempty(iroi)
                fprintf('WARNING: Failed to find the ROI %s in the %d-mm data of subject %s\n', ...
                        t_roi, xmm, sID);
                continue;
            end

            if isequal(meas, 'FA')
                fa.(grp)(h1, i2) = meanfa(iroi);
            elseif isequal(meas, 'L1')
                fa.(grp)(h1, i2) = meanL1(iroi);
            elseif isequal(meas, 'RD')
                fa.(grp)(h1, i2) = meanRD(iroi);
            else
                error('Unrecognized measure type: %s', meas);
            end
        end
    end
end
        
%% ROI-by-ROI comparisons
p_vals_byROI = nan(nROIs, 1);
t_vals_byROI = nan(nROIs, 1);
sig_vals_byROI = nan(nROIs, 1);

for i1 = 1 : nROIs
    xs_PWS = fa.PWS(i1, :);
    xs_PFS = fa.PFS(i1, :);
    
    [h, p, ci, stats] = ttest2(xs_PWS, xs_PFS);
    
    p_vals_byROI(i1) = p;
    t_vals_byROI(i1) = stats.tstat;
    sig_vals_byROI(i1) = -log10(p) * sign(t_vals_byROI(i1));
end


%% Print the significant differences
idx_sig = find(p_vals_byROI < TTEST_P_THRESH_UC);

for i1 = 1 : numel(idx_sig)
    fprintf('%s: p = %f; t = %f ', ...
            roi_names{idx_sig(i1)}, ...
            p_vals_byROI(idx_sig(i1)), ...
            t_vals_byROI(idx_sig(i1)))
        
    if t_vals_byROI(idx_sig(i1)) < 0
        fprintf('(PWS < PFS)\n');
    else
        fprintf('(PWS > PFS)\n');
    end
end

%% Prepare for behavioral correlation
roi_mean_fa = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    roi_mean_fa.(grp) = nan(numel(sIDs.(grp)), numel(ROIS_BEHAVCORR));
    
    for i2 = 1 : numel(ROIS_BEHAVCORR)
        t_roi = ROIS_BEHAVCORR{i2};
        
        i_roi = strmatch(t_roi, roi_names, 'exact');
        i_roi = fsic(roi_names, t_roi);
        
        if length(i_roi) == 1
            roi_mean_fa.(grp)(:, i2) = fa.(grp)(i_roi, :)';
        else
            fprintf(2, 'WARNING: Cannot find data for ROI %s for behavioral correlation.\n', ...
                    t_roi);
        end
    end
end
roi_mean_fa_2g = [roi_mean_fa.PWS; roi_mean_fa.PFS];

%% Correlation between FA and SSI4
[p_SSI4_corr_PWS, r_SSI4_corr_PWS] = corr_rois(roi_mean_fa.PWS, SSI4.PWS);

hashROIs.lh = {};
hashROIs.rh = {};
hashROISigns.lh = [];
hashROISigns.rh = [];

fprintf(1, '=== Significant correlations between ROI-mean FA and SSI4 ===\n');
for i1 = 1 : numel(ROIS_BEHAVCORR)
    if (p_SSI4_corr_PWS(i1) < LINCORR_P_THRESH_UC)
        fprintf(1, '%s: p = %.4f; r = %.3f\n', ...
                ROIS_BEHAVCORR{i1}, p_SSI4_corr_PWS(i1), r_SSI4_corr_PWS(i1));
            

        t_hemi = ROIS_BEHAVCORR{i1}(1 : 2);
        hashROIs.(t_hemi){end + 1} = ROIS_BEHAVCORR{i1}(4 : end);
        
        hashROISigns.(t_hemi)(end + 1) = ...
            (r_SSI4_corr_PWS(i1) > 0) * 2 - 1;
    end
end

%% Draw the aparc12 (SLaparc) ROI figures
% Create green-red colormap
cm0 = create_green_red_colormap();
cm = nan(size(cm0, 1) + 1, size(cm0, 2));
for i1 = 1 : size(cm, 2)
    cm(:, i1) = interp1(1 : size(cm0, 1), cm0(:, i1), ...
                        linspace(1, size(cm0, 1), size(cm, 1)));
end


for i1 = 1 : numel(hemis)
    hemi = hemis{i1};
    
    t_fillROIs = {};
    t_sigVals = [];
    for i2 = 1 : numel(roi_names)
        if isequal(roi_names{i2}(1 : 2), hemi)
            t_fillROIs{end + 1} = roi_names{i2}(4 : end);
            t_sigVals(end + 1) = sig_vals_byROI(i2);
        end
    end
    
    % -- Figure out the mapping between sig values and color values -- %
    max_abs_sig = max(abs(t_sigVals));
    t_fillClrs = cell(size(t_fillROIs));
    
    for i2 = 1 : numel(t_fillROIs)
        tpos = (t_sigVals(i2) + max_abs_sig) / 2 / max_abs_sig;
        tpos = round(1 + tpos * (size(cm, 1) - 1));
        t_fillClrs{i2} = cm(tpos, :);
    end
    
    % -- Determine which ROI names to highlight -- %
    hiliteROIs = t_fillROIs(find(abs(t_sigVals) > -log10(0.05)));
    hiliteROIClrs = repmat({hiliteROIClr}, 1, length(hiliteROIs));
    
    drawBndClrs = repmat({drawBndClr}, 1, length(hiliteROIs));
    
    % -- Carry out the drawing -- %
    imBase = imread(roiFigs.(hemi));
    imDiv = imread(roiFigsDiv.(hemi));
    imText = imread(roiFigsText.(hemi));
    
    imFASigMap.(hemi) = ...
        fill_aparc12_roi_fig(imBase, imDiv, imText, hemi, ...
                             t_fillROIs, 'fillClrs', t_fillClrs, ...
                             'clrNameROIs', hiliteROIs, hiliteROIClrs, ...
                             'drawBndROIs', hiliteROIs, drawBndClrs, ...
                             'hashROIs', hashROIs.(hemi), hashROISigns.(hemi))
                         
	% -- Draw the color bar with ticks -- %
    tickSigVals = 0 : 0.5 : max_abs_sig;
    tickSigVals = [-fliplr(tickSigVals(2 : end)), tickSigVals];
    tickSigIdx = (tickSigVals + max_abs_sig) / 2 / max_abs_sig;
    tickSigIdx = round(1 + tickSigIdx * (size(cm, 1) - 1));
    
    clrBarX = 400;
    clrBarY = 560;
    clrBarLenFact = 3;
    clrBarWidth = 20;
    
    imFASigMap.(hemi) = draw_colorbar(imFASigMap.(hemi), cm, ...
        [clrBarX, clrBarY], clrBarLenFact, clrBarWidth, ...
        'tickIdx', tickSigIdx);
    
	% -- Visualize -- %
	hf_faSigMap.(hemi) = figure('Name', sprintf('FA sig value map: %s', hemi), ...
                                'Color', defaultFigClr);
    
    imshow(imFASigMap.(hemi));
    hold on;
    
%     for j1 = 1 : length(tickSigIdx)
%         text(clrBarX + clrBarWidth + 10, ...
%             clrBarY + tickSigIdx(j1) * (clrBarLenFact + 1) / 2, ...
%             'a');
%     end
    
    % -- Save to image file -- %
    figFN = fullfile(figSaveDir, ...
                     sprintf('aparc12_%s.%dmm.%s.tif', meas, xmm, hemi));
    saveas(hf_faSigMap.(hemi) , figFN, 'tif');
    check_file(figFN);
    fprintf(1, 'Saved FA to image file: %s\n', figFN);
end




%% Perform correlations with behavioral measures
[p_rnSV_corr, r_rnSV_corr] = corr_rois(roi_mean_fa_2g, rnSV);
[p_rnSV_corr_PWS, r_rnSV_corr_PWS] = corr_rois(roi_mean_fa.PWS, rnSV.PWS);
[p_rnSV_corr_PFS, r_rnSV_corr_PFS] = corr_rois(roi_mean_fa.PFS, rnSV.PFS);

[p_hnSV_corr, r_hnSV_corr] = corr_rois(roi_mean_fa_2g, hnSV);
[p_hnSV_corr_PWS, r_hnSV_corr_PWS] = corr_rois(roi_mean_fa.PWS, hnSV.PWS);
[p_hnSV_corr_PFS, r_hnSV_corr_PFS] = corr_rois(roi_mean_fa.PFS, hnSV.PFS);

plot_corr(p_rnSV_corr, r_rnSV_corr, ...
          p_rnSV_corr_PWS, r_rnSV_corr_PWS, ...
          p_rnSV_corr_PFS, r_rnSV_corr_PFS, ...
          ROIS_BEHAVCORR, 'rnSV corr.', ...
          roi_mean_fa_2g, rnSV);
      
plot_corr(p_hnSV_corr, r_hnSV_corr, ...
          p_hnSV_corr_PWS, r_hnSV_corr_PWS, ...
          p_hnSV_corr_PFS, r_hnSV_corr_PFS, ...
          ROIS_BEHAVCORR, 'hnSV corr.', ...
          roi_mean_fa_2g, hnSV);
      


%% Comparison wih SfN 2011 poster

figure('Position', [50, 100, 800, 300], 'Color', 'w');
yLim = [0.25, 0.4];

fontSize = 15;
barW = 0.4;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
xTickLabel = {};
    
p_uc = nan(1, numel(ROIS_sel)); % Uncorrected p-value
for i1 = 1 : numel(ROIS_sel)
    t_ROI = ROIS_sel{i1};
    idx_roi = fsic(roi_names, t_ROI);
    t_FA_PWS = fa.PWS(idx_roi, :);
    t_FA_PFS = fa.PFS(idx_roi, :);
    bar(i1 - barW / 2, mean(t_FA_PFS), barW, 'FaceColor', 'none', ...
        'EdgeColor', colors.PFS, 'LineWidth', 1);
    hold on;
    plot(repmat(i1 - barW / 2, 1, 2), mean(t_FA_PFS) + [-1, 1] * ste(t_FA_PFS), ...
         'Color', colors.PFS, 'LineWidth', 1);
    bar(i1 + barW / 2, mean(t_FA_PWS), barW, 'FaceColor', 'none', ...
        'EdgeColor', colors.PWS, 'LineWidth', 1);
    plot(repmat(i1 + barW / 2, 1, 2), mean(t_FA_PWS) + [-1, 1] * ste(t_FA_PWS), ...
         'Color', colors.PWS, 'LineWidth', 1);

    xTickLabel{end + 1} = strrep(strrep(t_ROI, 'lh.', ''), 'rh.', '');

    [h_t, p_t] = ttest2(t_FA_PFS, t_FA_PWS);
    if p_t < 0.05
        plot(i1 + [-0.5, 0.5] * barW, repmat(yLim(2) - 0.075 * range(yLim), 1, 2), ...
             'b-', 'LineWidth', 1);
        plot(repmat(i1 - 0.5 * barW, 1, 2), yLim(2) - [0.1, 0.075] * range(yLim), ...
             'b-', 'LineWidth', 1);
        plot(repmat(i1 + 0.5 * barW, 1, 2), yLim(2) - [0.1, 0.075] * range(yLim), ...
             'b-', 'LineWidth', 1);
        plot(i1, yLim(2) - 0.025 * range(yLim), 'b*', 'MarkerSize', 10);        
    end
    p_uc(i1) = p_t; % Uncorrected p-value

    if i1 == 1
        text(i1 - barW, yLim(1) + 0.02 * range(yLim), 'PFS', ...
             'FontSize', fontSize, 'Color', colors.PFS);
        text(i1 + barW, yLim(1) + 0.02 * range(yLim), 'PWS', ...
             'FontSize', fontSize, 'Color', colors.PWS);
    end
end
set(gca, 'XLim', [0.25, numel(ROIS_sel) + 0.75]);
set(gca, 'XTick', 1 : numel(ROIS_sel), 'XTickLabel', xTickLabel);
set(gca, 'YLim', yLim);
ylabel('ROI FA (mean\pm1 SE)', 'FontSize', fontSize);      

%% Analysis and visualization
% meta_comp2(fa, sprintf('Mean FA'));

return        