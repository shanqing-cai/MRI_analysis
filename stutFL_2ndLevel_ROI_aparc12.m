function stutFL_2ndLevel_ROI_aparc12(mode_2l, p_thresh_2L_uc, bRandPerm)
%% Subject IDs
sIDs_all = {'S03', 'S06', 'S08', 'S09', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', ...
            'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S25', 'S26'  'S27', ...
            'S28', 'S29', 'S30', 'S31', 'S32', 'S33', 'S34', 'S35', 'S36', 'S37', ...
            'S38', 'S39'}; % S24 removed
        
%% Config: 2nd-level significance 
% p_thresh_2L_uc = 0.05; % uncorrected
FDR_q_thresh_2L = 0.05; % FDR
FWE_p_thresh_2L = 0.05; % FWE (rand perm)

p_thresh_bgc_uc = 0.05;  % uncorrected: between-group comparison
FDR_q_thresh_bgc = 0.05; % FDR: between-group comparison

FWE_N_PERM = 2e3;
    
%% Directories
roiDataDir = '/users/cais/STUT/analysis/aparc12_FL';

baseFigDir = '/users/cais/STUT/figures';

%% Group identities
sIDs.PWS = {};
sIDs.PFS = {};

for i1 = 1 : numel(sIDs_all)
    t_id = LUT_MRI_subjIDs(sIDs_all{i1}, 'inv');
    if isequal(t_id(1 : 3), 'PWS')
        sIDs.PWS{end + 1} = sIDs_all{i1};
    elseif isequal(t_id(1 : 3), 'PFS')
        sIDs.PFS{end + 1} = sIDs_all{i1};
    else
        error('Unrecognized subject ID: %s', t_id);
    end
end

%% Get all available ROIs
roiNames = {};
nROIs = 0;
groups = fields(sIDs);

for i1 = 1 : numel(groups)
    grp = groups{i1};
    
    for i2 = 1 : numel(sIDs.(grp))
        sID = sIDs.(grp){i2};
        
        dfs = dir(fullfile(roiDataDir, sID, 'roi_timecourses.*.mat'));
        
        if length(dfs) == 0
            error('WARNING: no roi_timecourses data files found for subject %s', sID);
        end
        
        for i3 = 1 : numel(dfs)
            dfn = fullfile(roiDataDir, sID, dfs(i3).name);
            load(dfn);  % gives roi_timecourses
            
            for k1 = 1 : roi_timecourses.nROIs
                t_roi = roi_timecourses.roiNames{k1};
                if isempty(fsic(roiNames, t_roi))
                    roiNames{end + 1} = t_roi;
                    nROIs = nROIs + 1;
                end
            end
        end          
    end
    
end

%% Cycle through all subjects and make tables of the t-values, for
%% second-level analysis
nSubjs.PFS = numel(sIDs.PFS);
nSubjs.PWS = numel(sIDs.PWS);

tVals = nan(nSubjs.PFS + nSubjs.PWS, nROIs);
conVals = nan(nSubjs.PFS + nSubjs.PWS, nROIs);  % contrasts

groups = {'PFS', 'PWS'};
sCnt = 1;
for i1 = 1 : numel(groups)
    grp = groups{i1};
    
    for i2 = 1 : numel(sIDs.(grp))
        sID = sIDs.(grp){i2};
        fprintf('Processing the data from subject %s\n', sID);
        
        dfs = dir(fullfile(roiDataDir, sID, 'roi_timecourses.*.mat'));
        
        if length(dfs) == 0
            fprintf('WARNING: no roi_timecourses data files found for subject %s', sID);
        end
        
        roitc_BL = []; % Baseline
        roitc_S = [];   % Speech
        
        for i3 = 1 : numel(dfs)
            dfn = fullfile(roiDataDir, sID, dfs(i3).name);
            fprintf('Loading data from: %s\n', dfn);
            load(dfn);  % gives roi_timecourses
            
            roitc_BL = [roitc_BL, roi_timecourses.dcv_data_BL];
            roitc_S = [roitc_S, roi_timecourses.dcv_data_S];
        end
        
        for k1 = 1 : roi_timecourses.nROIs
            t_data_BL = roitc_BL(k1, :);
            t_data_S = roitc_S(k1, :);

            [h, p, ci, stats] = ttest2(t_data_S, t_data_BL);
            
            t_roi = roi_timecourses.roiNames{k1};
            if isempty(fsic(roiNames, t_roi))
                error('Unexpected error: ROI = %s', t_ROI);
            else
                t_idx = fsic(roiNames, t_roi);
                tVals(sCnt, t_idx) = stats.tstat;
                conVals(sCnt, t_idx) = mean(t_data_S) - mean(t_data_BL);
            end
        end
        
        sCnt = sCnt + 1;
    end
    

end

%% Statistical analysis 1. Between-group comparison
% -- Uncorrected comparisons --
pVals_bgc_uc = nan(1, nROIs);
diffROIs_bgc_uc = {};
diffStr_bgc_uc = {};

for i1 = 1 : nROIs
    cons_PFS = conVals(1 : nSubjs.PFS, i1);
    cons_PWS = conVals(nSubjs.PFS + 1 : end, i1);
    
    [h, pVals_bgc_uc(i1)] = ttest2(cons_PWS, cons_PFS);
    if mean(cons_PWS) >  mean(cons_PFS)
        diffStr = 'PWS > PFS';
    else
        diffStr = 'PWS < PFS';
    end
    
    if pVals_bgc_uc(i1) < p_thresh_bgc_uc
        fprintf('Significant between-group diff. (p < %f): %s (p = %f, %s)\n', ...
                p_thresh_bgc_uc, roiNames{i1}, pVals_bgc_uc(i1), diffStr)
        diffROIs_bgc_uc{end + 1} = roiNames{i1};
        diffStr_bgc_uc{end + 1} = diffStr;
    end
end
fprintf('\n');

% -- FDR correction --
addpath('/users/cais/dtitools');
FDR_qVals_bgc = fdr1(pVals_bgc_uc);
for i1 = 1 : nROIs
    if FDR_qVals_bgc(i1) < FDR_q_thresh_bgc
        fprintf('Significant between-group diff (q < %f): %s (q = %f)\n', ...
                FDR_q_thresh_bgc, roiNames{i1}, FDR_qVals_bgc(i1), diffStr_bgc_uc{i1});
    end
end


%% Statistical analysis 2. Within-group 2L activation
activROIs_uc = struct;
activROIs_fdr = struct;
activROIs_fwe = struct;

pVals_2L_uc = struct;
FDR_qVals_2L = struct;
FWE_pVals_2L = struct;

nSig_2L_uc = struct;
nSig_2L_fdr = struct;
nSig_2L_fwe = struct;

for i1 = 1 : numel(groups) + 1
    if i1 <= numel(groups)
        grp = groups{i1}; 
    else
        grp = 'both';
    end
    
    activROIs_uc.(grp) = {};
    activROIs_fdr.(grp) = {};
    activROIs_fwe.(grp) = {};
    
    pVals_2L_uc.(grp) = nan(1, nROIs);
    FDR_qVals_2L.(grp) = nan(1, nROIs);
    FWE_pVals_2L.(grp) = nan(1, nROIs);
    
    nSig_2L_uc.(grp) = 0;
    nSig_2L_fdr.(grp) = 0;
    nSig_2L_fwe.(grp) = 0;
    
    if isequal(mode_2l, 't')
        if i1 == 1 % PFS        
            gtVals = tVals(1 : nSubjs.PFS, :);
        elseif i1 == 2
            gtVals = tVals(nSubjs.PFS + 1 : end, :);
        else
            gtVals = tVals;
        end
    elseif isequal(mode_2l, 'con')
        if i1 == 1 % PFS        
            gtVals = conVals(1 : nSubjs.PFS, :);
        elseif i1 == 2
            gtVals = conVals(nSubjs.PFS + 1 : end, :);
        else
            gtVals = conVals;
        end
    end
    
    % -- Uncorrected --
    for i2 = 1 : nROIs
        t_tVals = gtVals(:, i2);
        t_tVals = t_tVals(~isnan(t_tVals));
        
        [h, t_p_2L] = ttest(t_tVals);
        pVals_2L_uc.(grp)(i2) = t_p_2L;
        
        if t_p_2L < p_thresh_2L_uc
            activROIs_uc.(grp){end + 1} = roiNames{i2};
            fprintf(1, '%s (unc.): ROI %s: 2L p = %f (mean = %f)\n', ...
                    grp, roiNames{i2}, t_p_2L, mean(t_tVals));
            nSig_2L_uc.(grp) = nSig_2L_uc.(grp) + 1;
        end
    end
    
    fprintf('\n%s: TOTAL significant ROIs (unc.) = %d (%d / %d = %.2f%%)\n\n', ...
            grp, nSig_2L_uc.(grp), nSig_2L_uc.(grp), ...
            nROIs, nSig_2L_uc.(grp) / nROIs * 1e2);
    
    % -- FDR correction --
    addpath('/users/cais/dtitools');
    FDR_qVals_2L.(grp) = fdr1(pVals_2L_uc.(grp));
    for i2 = 1 : nROIs
        if FDR_qVals_2L.(grp)(i2) < FDR_q_thresh_2L
            activROIs_fdr.(grp){end + 1} = roiNames{i2};
            fprintf('%s (FDR): ROI %s: 2L q = %f\n', ...
                    grp, roiNames{i2}, FDR_qVals_2L.(grp)(i2));
            nSig_2L_fdr.(grp) = nSig_2L_fdr.(grp) + 1;
        end
    end
    
    fprintf('\n%s: TOTAL significant ROIs (FDR, q = %f) = %d (%d / %d = %.2f%%)\n\n', ...
            grp, FDR_q_thresh_2L, nSig_2L_fdr.(grp), nSig_2L_fdr.(grp), ...
            nROIs, nSig_2L_fdr.(grp) / nROIs * 1e2);
            
    % -- Random permutation FWE --
    if bRandPerm
        fprintf('Running random-permutation FWE on group %s...\n', grp);
        FWE_pVals_2L.(grp) = randperm_ttest_1samp(gtVals, FWE_N_PERM);
        for i2 = 1 : nROIs
            if FWE_pVals_2L.(grp)(i2) < FWE_p_thresh_2L
                activROIs_fwe.(grp){end + 1} = roiNames{i2};
                fprintf('%s (FWE, nperm = %d): ROI %s: 2L p = %f\n', ...
                        grp, FWE_N_PERM, ...
                        roiNames{i2}, FWE_pVals_2L.(grp)(i2));
                nSig_2L_fwe.(grp) = nSig_2L_fwe.(grp) + 1;
            end
        end

        fprintf('\n%s: TOTAL significant ROIs (FWE, nperm = %d, alpha = %f) = %d (%d / %d = %.2f%%)\n\n', ...
                grp, FWE_N_PERM, FWE_p_thresh_2L, ...
                nSig_2L_fwe.(grp), nSig_2L_fwe.(grp), ...
                nROIs, nSig_2L_fwe.(grp) / nROIs * 1e2);
    end
end

%% Write ROI sets to .mat file:
roiSet = struct;
dsFN = sprintf('/users/cais/STUT/scripts/activROIs_uc_thr%.3f.mat', p_thresh_2L_uc);
hemis = {'lh', 'rh'};

for i1 = 1 : numel(groups)
    grp = groups{i1};
    
    
    roiSet.(grp) = struct;
    for i2 = 1 : numel(hemis)
        hemi = hemis{i2};
        
        roiSet.(grp).(hemi) = {};
        
        for i3 = 1 : numel(activROIs_uc.(grp))
            if isequal(activROIs_uc.(grp){i3}(1 : 2), hemi)
                roiSet.(grp).(hemi){end + 1} = activROIs_uc.(grp){i3}(4 : end);                
            end
        end
    end
end

save(dsFN, 'roiSet');
if isfile(dsFN)
    fprintf(1, 'INFO: Saved ROI set to file %s\n', dsFN);
else
    error('Failed to write ROI set to file %s\n', dsFN);
end

%% Draw ROI figures
for i1 = 1 : numel(hemis)
    hemi = hemis{i1};
    
    roiFigs.(hemi) = fullfile(baseFigDir, sprintf('rois_%s_flat_SLaparc_noText.tif', hemi));
    roiFigsDiv.(hemi) = fullfile(baseFigDir, sprintf('rois_%s_flat_SLaparc_noText_div.tif', hemi));
    roiFigsText.(hemi) = fullfile(baseFigDir, sprintf('rois_%s_flat_SLaparc.tif', hemi));

    check_file(roiFigs.(hemi));
    check_file(roiFigsDiv.(hemi));
    check_file(roiFigsText.(hemi));
    
    im.(hemi) = imread(roiFigs.(hemi));
    imDiv.(hemi) = imread(roiFigsDiv.(hemi));
    imText.(hemi) = imread(roiFigsText.(hemi));
end

hf = [];
for i1 = 1 : numel(groups)
    grp = groups{i1};
    
    fillROIs = struct();
    for i2 = 1 : numel(hemis) + 1
        if i2 < 3
            hemi = hemis{i2};
            fillROIs.(hemi) = {};
            for i2 = 1 : numel(activROIs_uc.(grp))
                if isequal(activROIs_uc.(grp){i2}(1 : 3), [hemi, '_'])
                    fillROIs.(hemi){end + 1} = activROIs_uc.(grp){i2}(4 : end);
                end
            end
            vhemi = hemi;
        else
            hemi = 'union';
            vhemi = 'lh';
            fillROIs.(hemi) = union(fillROIs.(hemis{1}), fillROIs.(hemis{2}));
        end

        imFilled.(hemi) = fill_aparc12_roi_fig(im.(vhemi), imDiv.(vhemi), imText.(vhemi), ...
                                        vhemi, fillROIs.(hemi));    

        figFN = fullfile(baseFigDir, ...
                         sprintf('SpeechNetwork_stutFL_%s_%s_thr%.3f.tif', ...
                                 grp, hemi, p_thresh_2L_uc));

        hf(end + 1) = figure('Color', 'w', 'Name', ...
            sprintf('Speech network: %s, %s (n=%d)', grp, hemi, numel(fillROIs.(hemi))));
        imshow(imFilled.(hemi));
        saveas(hf(end), figFN, 'tif');
        check_file(figFN);
        fprintf('Saved ROI image from group %s, hemi %s to %s\n', ...
                grp, hemi, figFN);
    end
end

%% Save to data set
dsfn = fullfile('/users/cais/STUT/scripts', ...
                [mfilename, sprintf('_thr%.3f_ds.mat', p_thresh_2L_uc)]);
save(dsfn, 'activROIs_uc', 'activROIs_fdr', 'activROIs_fwe', ...
            'pVals_2L_uc', 'FDR_qVals_2L', 'FWE_pVals_2L', ...
            'nSig_2L_uc', 'nSig_2L_fdr', ...
            'p_thresh_2L_uc', 'FDR_q_thresh_2L', 'FWE_p_thresh_2L', ...
            'FWE_N_PERM', ...
            'sIDs');
        
if bRandPerm
    save(dsfn, 'nSig_2L_fwe', '-append');
end
fprintf('\nResults saved to mat file: %s\n', dsfn);
return

