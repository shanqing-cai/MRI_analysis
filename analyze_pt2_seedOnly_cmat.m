function analyze_pt2_seedOnly_cmat(hemi, netwName, meas, cmat_thresh, ...
                                   bFold, bReload, varargin)
%%
% Input: meas - measure to analyze {'tmn', 'wtn'}
%        cmat_thresh - threshold for binary global efficiency analysis
%        bFold - whether each individual subject's matrix should be
%                transposed and average with itself.
%        bReload - whether the time-consuming data reloading should be
%                  carried out
%
%%
sIDs.PWS = {'S01', 'S04', 'S06', 'S07', 'S08', 'S09', 'S10', 'S12', 'S15',  ...
            'S16', 'S17', 'S20', 'S21', 'S26', 'S28', 'S29', 'S33', 'S34', ...
            'S36', 'S37'};
            % Included S16 (But S16 should probably be kept in the group FA comparison) 
            % What's wrong with S21? Why was he excluded before, in
            % aparcSL_FA_analysis?
sIDs.PFS = {'S02', 'S03', 'S05', 'S11', 'S13', 'S14', 'S18', 'S19', 'S22', ...
            'S23', 'S25', 'S27', 'S30', 'S31', 'S32', 'S35', ...
            'S39'};
            % Left out S24, S38
            
TRACT_RES_DIR = '/users/cais/STUT/analysis/aparc12_tracts_pt2';
RSFC_BASE_DIR = '/users/cais/STUT/analysis';
RSFC_FILE_WC = 'corr_z_roi_aparc12.znbp2.noS24S38.%s.%s.mat';

p_thresh_unc = 0.01;
p_thresh_crl_unc = 0.005;

figSaveDir = '/users/cais/STUT/figures';

roiFigs.lh = '/users/cais/STUT/figures/rois_lh_flat_SLaparc.tif';
roiFigs.rh = '/users/cais/STUT/figures/rois_rh_flat_SLaparc.tif';
            
%%
bMaleOnly = ~isempty(fsic(varargin, 'maleOnly'));

%% Get speech-network ROI set
VALID_NETW_NAMES = {'speech', 'speech_PFS_lh', 'speec_PFS_rh', ...
                    'speech_PFS_lh', 'speech_PWS_rh', ...
                    'speech_2g_lh', 'speech_2g_rh'};

if length(netwName) > 6
    idxus = strfind(netwName, '_');
    if length(idxus) ~= 3
        error('Cannot find exactly three underlines in netwName: %s', netwName);
    end
    netwName0 = netwName(1 : idxus(3) - 1);
    p_thr = str2double(netwName(idxus(3) + 1 : end));
    
    if isempty(fsic(VALID_NETW_NAMES, netwName0))
        error('Unrecognized network name: %s', netwName);
    end
end

if ~isequal(netwName, 'speech')
    if ~isequal(netwName0(end - 1 : end), hemi)
        error('Hemisphere mismatch between netwName and hemi');
    end    
end
dsFN = sprintf('analyze_pt2_seedOnly_cmat_ds.%s.%s.mat', ...
                    netwName, hemi);

sprois = get_aparc12_cortical_rois(netwName, hemi);
nrois = length(sprois);

grps = fields(sIDs);

%% Optional: Load resting-state data 
bRSFC = ~isempty(fsic(varargin, '--rsfc'));
if bRSFC
    rsfcFN = fullfile(RSFC_BASE_DIR, sprintf(RSFC_FILE_WC, netwName, hemi));
    check_file(rsfcFN);
    load(rsfcFN); % gives corrDat
    
    if ~isequal(sprois, corrDat.spROIs)
        error('rsfc mode: mismatch between sprois and corrDat.spROIs');
    end
    
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        corrDat.mean_cm.(grp) = mean(corrDat.cms.(grp), 3);
    end
    
    % --- Take intersection of the subjects in the DTI and rsfc data sets --- %    
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        corrDat.sIDs.(grp) = cellstr(corrDat.sIDs.(grp));
        
        sIDs.(grp) = intersect(sIDs.(grp), corrDat.sIDs.(grp));
        if length(sIDs.(grp)) > length(corrDat.sIDs.(grp))
            error('There are subjects who are present in sIDs, but not in corrDat.sIDs');
        end
    end
end

%%
if bReload
    a_cmat = struct;    

    for i1 = 1 : numel(grps)
        grp = grps{i1};

        a_cmat.(grp) = nan(nrois, nrois, numel(sIDs.(grp)));

        for i2 = 1 : numel(sIDs.(grp))
            sID = sIDs.(grp){i2};
            fprintf(1, 'Loading data from subject (%s) %s...\n', grp, sID);

            mat_fn = fullfile(TRACT_RES_DIR, sID, ...
                        sprintf('connmats.%s.pt2.%s.mat', netwName, hemi));
            sdat = load(mat_fn);
            if isequal(meas, 'tmn')
                t_cmat = sdat.connmat_mean_norm;
            end
            
            if size(t_cmat, 1) ~= nrois
                error('Mismatch between nrois = %d and matrix size in: %s', ...
                      nrois, mat_fn);
            end
            
            a_cmat.(grp)(:, :, i2) = t_cmat;
    %         figplot(t_cmat_tmn(:), t_cmat_wtn(:), 'bo');
        end
    end

    save(dsFN, 'a_cmat');
else
    load(dsFN);
end

if bRSFC
    show_grp = 'PWS';
    show_sidx = 2;
    
    rs_tract_rsfc_corr = struct;
    ps_tract_rsfc_corr = struct;
    
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        
        rs_tract_rsfc_corr.(grp) = nan(1, numel(sIDs.(grp)));
        ps_tract_rsfc_corr.(grp) = nan(1, numel(sIDs.(grp)));
        
        for i2 = 1 : numel(sIDs.(grp))
            x = a_cmat.(grp)(:, :, i2);
            x = x(:);
            y = corrDat.cms.(grp)(:, :, i2);
            y = y(:);
            [~, r_tract_rsfc_corr, p_tract_rsfc_corr] = lincorr(x, y);
            
            rs_tract_rsfc_corr.(grp)(i2) = r_tract_rsfc_corr;
            ps_tract_rsfc_corr.(grp)(i2) = p_tract_rsfc_corr;

            if isequal(grp, show_grp) && i2 == show_sidx;  
                figure('Color', 'w');
                semilogx(x, y, 'o');
                xlabel('Normalized tract density');
                ylabel('Fisher-transformed BOLD correlation coefficient ');
                xs = get(gca, 'XLim'); 
                ys = get(gca, 'YLim');
                text(xs(1) + 0.0001 * range(xs), ys(2) - 0.06 * range(ys), ...
                     sprintf('R = %.3f; p = %.3f', r_tract_rsfc_corr, p_tract_rsfc_corr));
            end
        end
    end
    
    figure('Color', 'w', 'Position', [200, 200, 400, 300]);
    errorbar([1, 2], [mean(rs_tract_rsfc_corr.PFS), mean(rs_tract_rsfc_corr.PWS)], ...
             [ste(rs_tract_rsfc_corr.PFS), ste(rs_tract_rsfc_corr.PWS)]);
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
    ylabel('Correlation coefficient (R) (Mean \pm 1 SEM)');
end

%% Exclude subjects based on gender (optional)
if bMaleOnly
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        
        bKeep = ones(1, numel(sIDs.(grp)));
        
        for i2 = 1 : numel(sIDs.(grp))
            t_gend = get_subj_gender(sIDs.(grp){i2});
            if isequal(t_gend, 'Female')
                bKeep(i2) = 0;
            end
        end
        
        sIDs.(grp) = sIDs.(grp)(find(bKeep));
        
        a_cmat.(grp) = a_cmat.(grp)(:, :, find(bKeep));        
    end
end

%% 
if bFold
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        
        a_cvec.(grp) = [];
        for i2 = 1 : size(a_cmat.(grp), 3)
            a_cmat.(grp)(:, :, i2) = ...
                (a_cmat.(grp)(:, :, i2) + a_cmat.(grp)(:, :, i2)') / 2;
            
            t_mat = triu(a_cmat.(grp)(:, :, i2));
            t_vec = [];
            for k1 = 1 : nrois
                t_vec = [t_vec; t_mat(1 : k1 - 1, k1)];
            end
            
            a_cvec.(grp) = [a_cvec.(grp); t_vec'];                     
        end
    end
        
    save(dsFN, 'a_cvec', '-append');
end

%% Create the seed-by-seed vectors
a_seedVec = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    for i2 = 1 : numel(sIDs.(grp))
        a_seedVec.(grp) = squeeze(nanmean(a_cmat.(grp), 2));
    end
end

%% Compare the seed-by-seed mean vectors
sv_p_rs = nan(1, numel(sprois));
sv_sgn_rs = nan(1, numel(sprois));
for i1 = 1 : numel(sprois)
    dat = a_seedVec;    
    
    for i2 = 1 : numel(sprois)
        t_vec_PWS = dat.PWS(i2, :);
        t_vec_PFS = dat.PFS(i2, :);
        
        [sv_p_rs(i2), ~, stats] = ranksum(t_vec_PFS(:), t_vec_PWS(:));
        if median(t_vec_PWS) < median(t_vec_PFS)
            sv_sgn_rs(i2) = -1;
        else
            sv_sgn_rs(i2) = 1;
        end
    end
end

%% Load the SSI4 scores
SSI4 = nan(1, length(sIDs.PWS));
for i1 = 1 : numel(sIDs.PWS)
    SSI4(i1) = ds_SSI4(sIDs.PWS{i1});
end

%% Load the EH_comp_300 scores 
%%     (magnitude of compensation to EH F1 perturbation, at 300 ms following pert. onset)
EH_comp = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    EH_comp.(grp) = nan(1, numel(sIDs.(grp)));
    
    for i2 = 1 : numel(sIDs.(grp))
        EH_comp.(grp)(i2) = get_qdec_measure(sIDs.(grp){i2}, 'EH_comp_300');
    end
end

%% Load the rnSV values
rnSV = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    rnSV.(grp) = nan(1, numel(sIDs.(grp)));
    
    for i2 = 1 : numel(sIDs.(grp))
        rnSV.(grp)(i2) = get_qdec_measure(sIDs.(grp){i2}, 'rnSV');
    end
end



%% Element-by-element between-group comparisons and correlations
% if isequal(meas, 'tmn')
%     a_cmat = a_cmat_tmn;
% elseif isequal(meas, 'wtn')
%     a_cmat = a_cmat_wtn;
% else
%     error('Unrecognized value of meas: %s', meas);
% end

% p_t = nan(nrois, nrois); % p-values from t-tests
p_rs = nan(nrois, nrois); % p-values from rank-sum tests
sgn_rs = nan(nrois, nrois);

[p_rs, sgn_rs, ...
    p_SSI4_spr, rho_SSI4_spr, ...
    p_EHcomp_spr, rho_EHcomp_spr, ...
    p_rnSV_spr, rho_rnSV_spr] = ...
        compare_cmats(a_cmat, bFold, 'ranksum', ...
                      SSI4, EH_comp, rnSV);
                  
if bRSFC
    [p_rs_rsfc, sgn_rs_rsfc, ...
        p_SSI4_spr_rsfc, rho_SSI4_spr_rsfc, ...
        p_EHcomp_spr_rsfc, rho_EHcomp_spr_rsfc, ...
        p_rnSV_spr_rsfc, rho_rnSV_spr_rsfc] = ...
            compare_cmats(corrDat.cms, bFold, 'ranksum', ...
                          SSI4, EH_comp, rnSV);
end

pnSigs(1) = numel(find(p_rs(:) < 0.05 & sgn_rs(:) == 1));
pnSigs(2) = numel(find(p_rs(:) < 0.05 & sgn_rs(:) == -1));
p_bn_bgd = myBinomTest(pnSigs(1), sum(pnSigs), 0.5, 'Two');
fprintf(1, 'Binomial test on the between-group difference results (two-tailed): p = %e\n', ...
        p_bn_bgd);
    
pnSigs(1) = numel(find(p_SSI4_spr(:) < 0.05 & rho_SSI4_spr(:) > 0));
pnSigs(2) = numel(find(p_SSI4_spr(:) < 0.05 & rho_SSI4_spr(:) < 0));
p_bn_corrSSI4 = myBinomTest(pnSigs(1), sum(pnSigs), 0.5, 'Two');
fprintf(1, 'Binomial test on the SSI-4 correlation results (two-tailed): p = %e\n', ...
        p_bn_corrSSI4);

%% --- Perform random permutation --- %% 
if ~isempty(fsic(varargin, '--randPerm'))
    nRandPerm = varargin{fsic(varargin, '--randPerm') + 1};
    pnSigs = nan(nRandPerm, 2); % Numbers of positive and negative differences
    
    % aa_cmat = cat(3, a_cmat.PFS, a_cmat.PWS);
    % nts = length(sIDs.PFS) + length(sIDs.PWS);
    rp_p_rs = nan(nrois, nrois, nRandPerm);
    rp_sgn_rs = nan(nrois, nrois, nRandPerm);
    
    bPermSSI4 = ~isempty(fsic(varargin, '--randPermSSI4'));
    
    pnSigs_SSI4 = nan(nRandPerm, 2);
    rp_p_SSI4_spr = nan(nrois, nrois, nRandPerm);
    rp_rho_SSI4_spr = nan(nrois, nrois, nRandPerm);
    
%     for i1 = 1 : nRandPerm
    parfor i1 = 1 : nRandPerm
        if mod(i1, 100) == 0
            fprintf(1, 'Performing random permutation %d of %d...\n', i1, nRandPerm);
        end
        
        if ~bPermSSI4
            [rp_p_rs(:, :, i1), rp_sgn_rs(:, :, i1)] = ...
                compare_cmats(a_cmat, bFold, 'ranksum', '--randPerm');
        else
            [rp_p_rs(:, :, i1), rp_sgn_rs(:, :, i1), ...
                rp_p_SSI4_spr(:, :, i1), rp_rho_SSI4_spr(:, :, i1), ...
                ~, ~, ~, ~] = ...
                compare_cmats(a_cmat, bFold, 'ranksum', ...
                              SSI4, EH_comp, rnSV, '--randPerm');
        end

%         pnSigs(i1, 1) = numel(find(rp_p_rs(:, :, i1) < 0.05 & rp_sgn_rs(:, :, i1) == 1));
%         pnSigs(i1, 2) = numel(find(rp_p_rs(:, :, i1) < 0.05 & rp_sgn_rs(:, :, i1) == -1));
        pnSigs(i1, :) = [numel(find(rp_p_rs(:, :, i1) < 0.05 & rp_sgn_rs(:, :, i1) == 1)), ...
                         numel(find(rp_p_rs(:, :, i1) < 0.05 & rp_sgn_rs(:, :, i1) == -1))];
        if bPermSSI4
            pnSigs_SSI4(i1, :) = [numel(find(rp_p_SSI4_spr(:, :, i1) < 0.05 & rp_rho_SSI4_spr(:, :, i1) > 0)), ...
                                  numel(find(rp_p_SSI4_spr(:, :, i1) < 0.05 & rp_rho_SSI4_spr(:, :, i1) < 0))];
        end
    end
    rp_pnRatios = (pnSigs(:, 1) + 1) ./ (pnSigs(:, 2) + 1);
    

    pnSigs = nan(1, 2);
    pnSigs(1) = numel(find(p_rs(:) < 0.05 & sgn_rs(:) == 1));
    pnSigs(2) = numel(find(p_rs(:) < 0.05 & sgn_rs(:) == -1));
    pnRatios = ((pnSigs(1) + 1) ./ (pnSigs(2) + 1));
    
    fprintf(1, 'Permutation test (N=%d) on the between-group difference results (one-tailed): p = %e\n', ...
        nRandPerm, length(find(rp_pnRatios <= pnRatios)) / nRandPerm);
   
    if bPermSSI4
        fprintf(1, 'Permutation test (N=%d) on the correlation with SSI4 (one-tailed):\n', nRandPerm);
        % --- Test the significant of number of significant differences ---
        rp_totSig_SSI4 = sum(pnSigs_SSI4, 2);
        totSig_SSI4 = numel(find(p_SSI4_spr(:) < 0.05));
        
        fprintf(1, '\tOn the number of significant correlations (regardless of sign): p = %e\n', ...
                length(find(rp_totSig_SSI4 > totSig_SSI4)) / nRandPerm);
        
        rp_pnRatios_SSI4 = (pnSigs_SSI4(:, 1) + 1) ./ (pnSigs_SSI4(:, 2) + 1);
        pnSigs_SSI4 = nan(1, 2);
        pnSigs_SSI4(1) = numel(find(p_SSI4_spr(:) < 0.05 & rho_SSI4_spr(:) > 0));
        pnSigs_SSI4(2) = numel(find(p_SSI4_spr(:) < 0.05 & rho_SSI4_spr(:) < 0));
        pnRatios_SSI4 = ((pnSigs_SSI4(1) + 1) ./ (pnSigs_SSI4(2) + 1));
        
        
        fprintf(1, '\tBias toward negative correlations: p = %e\n', ...
                length(find(rp_pnRatios_SSI4 <= pnRatios_SSI4)) / nRandPerm);
    end
end

if ~isempty(fsic(varargin, 'randPermPWS'))
    nRandPermPWS = varargin{fsic(varargin, 'randPermPWS') + 1};
    pnSigs = nan(nRandPermPWS, 2); % Numbers of positive and negative differences
    % aa_cmat = cat(3, a_cmat.PFS, a_cmat.PWS);
    % nts = length(sIDs.PFS) + length(sIDs.PWS);
    
    rp_p_SSI4_spr = nan(nrois, nrois, nRandPermPWS);
    rp_rho_SSI4_spr = nan(nrois, nrois, nRandPermPWS);
    
    parfor i1 = 1 : nRandPermPWS
        if mod(i1, 10) == 0
            fprintf(1, 'Performing random permutation (PWS) %d of %d...\n', i1, nRandPermPWS);
        end
        [~, ~, rp_p_SSI4_spr(:, :, i1), rp_rho_SSI4_spr(:, :, i1), ~, ~, ~, ~] = ...
            compare_cmats(a_cmat, bFold, 'ranksum', SSI4, EH_comp, rnSV, '--randPermPWS');

%         pnSigs(i1, 1) = numel(find(rp_p_SSI4_spr(:, :, i1) < 0.05 & rp_rho_SSI4_spr(:, :, i1) > 0));
%         pnSigs(i1, 2) = numel(find(rp_p_SSI4_spr(:, :, i1) < 0.05 & rp_rho_SSI4_spr(:, :, i1) < 0));
        pnSigs(i1, :) = [numel(find(rp_p_SSI4_spr(:, :, i1) < 0.05 & rp_rho_SSI4_spr(:, :, i1) > 0)), ...
                         numel(find(rp_p_SSI4_spr(:, :, i1) < 0.05 & rp_rho_SSI4_spr(:, :, i1) < 0))];
    end
    rp_pnRatios = (pnSigs(:, 1) + 1) ./ (pnSigs(:, 2) + 1);

    pnSigs = nan(1, 2);
    pnSigs(1) = numel(find(p_SSI4_spr(:) < 0.05 & rho_SSI4_spr(:) > 0));
    pnSigs(2) = numel(find(p_SSI4_spr(:) < 0.05 & rho_SSI4_spr(:) < 0));
    pnRatios = ((pnSigs(1) + 1) ./ (pnSigs(2) + 1));
    
    fprintf(1, 'Permutation test (N=%d) on the SSI-4 correlation results (one-tailed): p = %e\n', ...
        nRandPermPWS, length(find(rp_pnRatios <= pnRatios)) / nRandPermPWS);
    end
% --- ~Perform random permutation --- % 

%% FDR
a_p_rs = p_rs(:);
a_p_rs = a_p_rs(~isnan(a_p_rs));

%% Print results from the between-group comparison
fprintf(1, '=== Significant differences from rank-sum test ===\n');
for i1 = 1 : nrois
    for i2 = 1 : nrois
        if p_rs(i1, i2) < p_thresh_unc
            med_PFS = median(squeeze(a_cmat.PFS(i1, i2, :)));
            med_PWS = median(squeeze(a_cmat.PWS(i1, i2, :)));
            
            if med_PWS < med_PFS
                diffDir = '<';
            else
                diffDir = '>';
            end
            
            % Look at SSI4 correlation
            if p_SSI4_spr(i1, i2) < 0.05
                if rho_SSI4_spr(i1, i2) < 0
                    crlStr = sprintf('-: %.3f', p_SSI4_spr(i1, i2));
                else
                    crlStr = sprintf('+: %.3f', p_SSI4_spr(i1, i2));                
                end
            else
                crlStr = '';
            end
            
            fprintf(1, '%s -> %s:\tmed = %f (PWS) %s %f (PFS); p = %e [%s]\n', ...
                    sprois{i1}, sprois{i2}, ...
                    med_PWS, diffDir, med_PFS, p_rs(i1,i2), crlStr);
        end
    end
end
fprintf(1, '\n');

%% Print results from the correlation with SSI4
fprintf(1, '=== Significant correlatoins with SSI4 (PWS only) ===\n');
for i1 = 1 : nrois
    for i2 = 1 : nrois
        if p_SSI4_spr(i1, i2) < p_thresh_crl_unc
            if rho_SSI4_spr(i1, i2) < 0
                corrDir = '-';
            else
                corrDir = '+';
            end
            
            fprintf(1, '%s -> %s: p = %e (%s)\n', ...
                    sprois{i1}, sprois{i2}, ...
                    p_SSI4_spr(i1,i2), corrDir);
        end
    end
end
fprintf(1, '\n');

%% BCT (Graph theory) analysis: Binary global efficiency
efb = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    efb.(grp) = nan(numel(sIDs.(grp)), 1);
    
    for i2 = 1 : numel(sIDs.(grp))
        binMat = eye(nrois);        
        binMat(a_cmat.(grp)(:, :, i2) > cmat_thresh) = 1;       
        
        efb.(grp)(i2) = efficiency_bin(binMat);
    end
end

p_bge = ranksum(efb.PWS, efb.PFS);

fprintf(1, '=== Binary global efficiency (cmat_thresh = %f) ===\n', cmat_thresh);
fprintf(1, 'Between-group: p = %e\n\n', p_bge);

%% BCT (Graph theory) analysis: Weighted global efficiency
efw = struct;

% weiMax_PWS = max(a_cmat.PWS(:));
% weiMax_PFS = max(a_cmat.PFS(:));
% weiMax = max([weiMax_PWS, weiMax_PFS]);
weiMax = 150;

for i1 = 1 : numel(grps)
    grp = grps{i1};
    efw.(grp) = nan(numel(sIDs.(grp)), 1);
    
    for i2 = 1 : numel(sIDs.(grp))
        weiMat = a_cmat.(grp)(:, :, i2) / weiMax;
        weiMat(weiMat > 1) = 1;
        
        efw.(grp)(i2) = efficiency_wei(weiMat);
    end
end

p_wge = ranksum(efw.PWS, efw.PFS);
[rho_spr_wge_SSI4, ~, p_spr_wge_SSI4] = spear(efw.PWS, SSI4(:));
[k_lc_wge_SSI4, r2_lc_wge_SSI4, p_lc_wge_SSI4] = lincorr(efw.PWS, SSI4(:));
r_lc_wge_SSI4 = sign(k_lc_wge_SSI4(2)) * sqrt(r2_lc_wge_SSI4);

fprintf(1, '=== Weighted global efficiency (weiMax = %f) ===\n', weiMax);
fprintf(1, 'Between-group: p = %e\n\n', p_wge);
fprintf(1, 'Spearman correl. with SSI4 (PWS): rho = %f, p = %e\n', ...
        rho_spr_wge_SSI4, p_spr_wge_SSI4);
fprintf(1, 'Linear correl. with SSI4 (PWS): r = %f, p = %e\n\n', ...
        r_lc_wge_SSI4, p_lc_wge_SSI4);   

%% Compute and visualize the average cmats and significance of differences

figSize = 800;
verticalPadding = 2.5;
horizontalPadding = 1.5;

mn_cmat = struct;


for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    mn_cmat.(grp) = mean(a_cmat.(grp), 3);
    
    
    figure('Position', [100, 300, figSize, figSize], 'Color', 'w', ...
           'Name', ['Connecivity matrix: ', grp]);
    axis tight;    
        
%     subplot(1, 2, i1);
    imagesc(mn_cmat.(grp));
    hold on;
    colormap bone;
    axis square;
    
%     title(grp);
    colorbar;
    
    set(gca, 'XTick', [], 'YTick', []);
    xs = get(gca, 'XLim');
    ys = get(gca, 'YLim');
    
    % --- Labels for columns --- %
    ht_cols = nan(1, numel(sprois));
    ht_cols_top = nan(1, numel(sprois));
    for k1 = 1 : numel(sprois)
        ht_cols(k1) = text(k1, numel(sprois) + verticalPadding, ...
                           strrep(sprois{k1}, [hemi, '_'], ''), ...
                           'FontSize', 12);
        ht_cols_top(k1) = text(k1, 0, ...
                           strrep(sprois{k1}, [hemi, '_'], ''), ...
                           'FontSize', 12);
        set(ht_cols(k1), 'rotation', 90);
        set(ht_cols_top(k1), 'rotation', 90);
    end
    
    % --- Labels for rows --- %
    ht_rows = nan(1, numel(sprois));
    for k1 = 1 : numel(sprois)
        ht_rows(k1) = text(-horizontalPadding, k1, ...
                           strrep(sprois{k1}, [hemi, '_'], ''), ...
                           'FontSize', 12);
    end
       
    plot(xs, ys, '-', 'Color', [0.5, 0.5, 0.5]);
end

%% Draw ROI figures
if bFold
    drawRoiPairs = struct;

    for i1 = 1 : numel(grps)
        grp = grps{i1};
        drawRoiPairs.(grp) = cell(0, 3);

        for k1 = 1 : numel(sprois)
            for k2 = 1 : numel(sprois)
                if bFold && (k2 > k1)
                    continue;
                end
                
                if mn_cmat.(grp)(k1, k2) > cmat_thresh / 5 % WARNING: ad hoc!
                    drawRoiPairs.(grp) = [drawRoiPairs.(grp); ...
                                          {sprois{k1}, sprois{k2}, mn_cmat.(grp)(k1, k2)}];
                end
            end
        end
        
        draw_aparc12_roi_arrows(roiFigs.(hemi), drawRoiPairs.(grp), ...
                                [20, 50, 100], ...
                                {[0, 0.5, 0], 'b', 'r'});
        title(grp);
        
        figFN = fullfile(figSaveDir, ...
            sprintf('aparc12_pt2_seedOnly_roiGraph.%s.%s.%s.%s.tif', ...
                    netwName, hemi, meas, grp));
        if bMaleOnly
            figFN = strrep(figFN, '.tif', '.maleOnly.tif');
        end
        saveas(gcf, figFN, 'tif');
        fprintf(1, 'Saved to image file: %s\n', figFN);
    end
end

%% --- Connectivity matrix difference --- %%
cellShift = 0.1;
visParams = [figSize, verticalPadding, horizontalPadding, cellShift];

sig_rs = sgn_rs .* -log10(p_rs);
sig_rs(isnan(sig_rs)) = 0;

[nSigs_bgd, sigConns_bgd] = ...
    show_2d_mat(sig_rs, sprois, hemi, ...
                'Connectivity matrix difference', visParams);
if bRSFC
    show_2d_mat(corrDat.ttest_sig, sprois, hemi, ...
                'rsFMRI connectivity matrix difference', visParams);
end

figFN = fullfile(figSaveDir, ...
    sprintf('aparc12_pt2_seedOnly_cmat_diff.%s.%s.%s.tif', ...
            netwName, hemi, meas));
if bMaleOnly
    figFN = strrep(figFN, '.tif', '.maleOnly.tif');
end
saveas(gcf, figFN, 'tif');
fprintf(1, 'Saved to image file: %s\n', figFN);

%% --- Connectivity matrix correlation with SSI4 --- %%
sig = sign(rho_SSI4_spr) .* -log10(p_SSI4_spr);
sig(isnan(sig)) = 0;

[nSigs_corrSSI4, sigConns_corrSSI4] = ...
    show_2d_mat(sig, sprois, hemi, ...
                'Connectivity correlation with SSI4', visParams);
% nSigs_corrSSI4 = numel(sigConns_corrSSI4);
disp(['Tractography connectivity correlation with SSI4: nSigs = ']);
disp(nSigs_corrSSI4);

if bRSFC
    sig_rsfc = sign(rho_SSI4_spr_rsfc) .* -log10(p_SSI4_spr_rsfc);
    sig_rsfc(isnan(sig_rsfc)) = 0.0;
    [nSigs_corrSSI4_rsfc, sigConns_corrSSI4_rsfc] = ...
        show_2d_mat(sig_rsfc, sprois, hemi, ...
                    'rsFMRI connectivity correlation with SSI4', visParams);
end

% --- Find the intersect of the sets of connections with significant bgd and
% those with significant correlation with SSI-4 scores ---
intersectConns = intersect_cells(sigConns_bgd, sigConns_corrSSI4);
fprintf(1, '%d ROIs with both significant between-group difference and significant correlation with SSI-4\n', ...
        length(intersectConns));
if length(intersectConns) > 0
    for i1 = 1 : length(intersectConns)
        fprintf(1, '\t%s <-> %s\n', intersectConns{i1}{1}, intersectConns{i1}{2});
    end
end


figFN = fullfile(figSaveDir, ...
    sprintf('aparc12_pt2_seedOnly_cmat_SSI4_spr.%s.%s.%s.tif', ...
            netwName, hemi, meas));
if bMaleOnly
    figFN = strrep(figFN, '.tif', '.maleOnly.tif');
end
saveas(gcf, figFN, 'tif');
fprintf(1, 'Saved to image file: %s\n', figFN);

%% -- Connectivity matrix correlation with EH_comp --- %
sig = sign(rho_EHcomp_spr) .* -log10(p_EHcomp_spr);
sig(isnan(sig)) = 0;

show_2d_mat(sig, sprois, hemi, ...
            'Connectivity correlation with EHcomp', visParams);
        
figFN = fullfile(figSaveDir, ...
    sprintf('aparc12_pt2_seedOnly_cmat_EHcomp_spr.%s.%s.%s.tif', ...
            netwName, hemi, meas));
if bMaleOnly
    figFN = strrep(figFN, '.tif', '.maleOnly.tif');
end
saveas(gcf, figFN, 'tif');
fprintf(1, 'Saved to image file: %s\n', figFN);

%% -- Connectivity matrix correlation with rnSV --- %
sig = sign(rho_rnSV_spr) .* -log10(p_rnSV_spr);
sig(isnan(sig)) = 0;

show_2d_mat(sig, sprois, hemi, ...
            'Connectivity correlation with rnSV', visParams);
        
figFN = fullfile(figSaveDir, ...
    sprintf('aparc12_pt2_seedOnly_cmat_rnSV_spr.%s.%s.%s.tif', ...
            netwName, hemi, meas));
if bMaleOnly
    figFN = strrep(figFN, '.tif', '.maleOnly.tif');
end
saveas(gcf, figFN, 'tif');
fprintf(1, 'Saved to image file: %s\n', figFN);

return