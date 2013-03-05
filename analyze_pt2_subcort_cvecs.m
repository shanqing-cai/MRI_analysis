function analyze_pt2_subcort_cvecs(hemi, scROI, segType, meas, ...
                                   cmat_thresh, bReload, varargin)
%% CONFIG: Subjects
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
            
%% CONFIG: Directory and file names
TRACT_RES_DIR = '/users/cais/STUT/analysis/aparc12_tracts_pt2';

figSaveDir = '/users/cais/STUT/figures';

dsFN = sprintf('analyze_pt2_subcort_cvecs_ds.%s.%s.seg%s.mat', hemi, scROI, segType);

SEG_TYPE_NUMSEGS.A = 7;
SEG_TYPE_NUMSEGS.A2_5 = 1;
SEG_TYPE_NUMSEGS.A2_4 = 1;

%% CONFIG: Statistical thresholds
p_thresh_unc = 0.05;
p_thresh_crl_unc = 0.005;

%%
bMaleOnly = ~isempty(fsic(varargin, 'maleOnly'));

%%
sprois = get_aparc12_cortical_rois('speech', hemi);
% sprois = get_aparc12_cortical_rois(hemi);
nrois = length(sprois);

grps = fields(sIDs);

a_cmat = struct;    
if length(strfind(segType, '-')) == 1
    segType1 = strrep(segType, '-', '_');
    nSegs = SEG_TYPE_NUMSEGS.(segType1);
else
    nSegs = SEG_TYPE_NUMSEGS.(segType);
end

if bReload
    for i1 = 1 : numel(grps)
        grp = grps{i1};

        a_cmat.(grp) = nan(nSegs, nrois, numel(sIDs.(grp)));

        for i2 = 1 : numel(sIDs.(grp))
            sID = sIDs.(grp){i2};
            fprintf(1, 'Loading data from subject (%s) %s...\n', grp, sID);

            mat_fn = fullfile(TRACT_RES_DIR, sID, ...
                        sprintf('connmats.pt2.%s.%s.seg%s.mat', ...
                                hemi, scROI, segType));
            sdat = load(mat_fn);
            
            if isequal(meas, 'tmn')
                t_cmat = sdat.connmat_mean_norm;
            else
                error('Unrecognized meas name: %s', meas);
            end
            
            % --- Find out the indices of the each of the sprois --- %
            idxs = nan(1, length(sprois));
            
            for k1 = 1 : length(sprois)
                idxs(k1) = strmatch(sprois{k1}, sdat.h_rois, 'exact');
            end
            
            a_cmat.(grp)(:, :, i2) = t_cmat(:, idxs);
                        
%             figplot(t_cmat_tmn(:), t_cmat_wtn(:), 'bo');
        end
    end

    save(dsFN, 'a_cmat');
    fprintf(1, 'INFO: saved data to file: %s\n', dsFN);
else
    load(dsFN);
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
        
        a_cmat_tmn.(grp) = a_cmat_tmn.(grp)(:, :, find(bKeep));
        a_cmat_wtn.(grp) = a_cmat_wtn.(grp)(:, :, find(bKeep));
    end
end

%% Load the SSI4 scores
SSI4 = nan(1, length(sIDs.PWS));
for i1 = 1 : numel(sIDs.PWS)
    SSI4(i1) = ds_SSI4(sIDs.PWS{i1});
end


%% Element-by-element between-group comparisons and correlations
p_t = nan(nSegs, nrois); % p-values from t-tests
p_rs = nan(nSegs, nrois); % p-values from rank-sum tests
sgn_rs = nan(nSegs, nrois);

rho_SSI4_spr = nan(nrois, nrois); 
p_SSI4_spr = nan(nrois, nrois); % p-values from Spearman correlation with SSI4 (PWS only)
% 
% rho_EHcomp_spr = nan(nrois, nrois);
% p_EHcomp_spr = nan(nrois, nrois);
% 
% rho_rnSV_spr = nan(nrois, nrois);
% p_rnSV_spr = nan(nrois, nrois);

for i1 = 1 : nSegs
    for i2 = 1 : nrois
        v_PFS = squeeze(a_cmat.PFS(i1, i2, :));
        v_PWS = squeeze(a_cmat.PWS(i1, i2, :));
        
        [foo, p_t(i1, i2)] = ttest2(v_PFS, v_PWS);
        [p_rs(i1, i2), foo, stats] = ranksum(v_PFS, v_PWS);
        if median(v_PWS) < median(v_PFS)
            sgn_rs(i1, i2) = -1;
        else
            sgn_rs(i1, i2) = 1;
        end
        
%         % --- Correlation with SSI4 --- %
        [rho_SSI4_spr(i1, i2), foo, p_SSI4_spr(i1, i2)] = ...
            spear(v_PWS, SSI4(:));
%         
%         % --- Correlation with EH_comp --- %
%         v_2grp = [v_PWS; v_PFS];
%         EHcomp_2grp = [EH_comp.PWS, EH_comp.PFS];
%         [rho_EHcomp_spr(i1, i2), ~, p_EHcomp_spr(i1, i2)] = ...
%             spear(v_2grp(:), EHcomp_2grp(:));
%         
%         % --- Correlation with rnSV --- %
%         rnSV_2grp = [rnSV.PWS, rnSV.PFS];
%         [rho_rnSV_spr(i1, i2), ~, p_rnSV_spr(i1, i2)] = ...
%             spear(v_2grp(:), rnSV_2grp(:));
    end
end

sig_rs = -log10(p_rs) .* sgn_rs;

%% Print results from the between-group comparison
fprintf(1, '=== Significant differences from rank-sum test (p_thresh_unc = %f) ===\n', ...
        p_thresh_unc);
for i1 = 1 : nSegs
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
            
            fprintf(1, 'Segment #%d -> %s:\tmed = %f (PWS) %s %f (PFS); p = %e [%s]\n', ...
                    i1, sprois{i2}, ...
                    med_PWS, diffDir, med_PFS, p_rs(i1,i2), crlStr);
%             fprintf(1, 'Segment #%d -> %s:\tmed = %f (PWS) %s %f (PFS); p = %e\n', ...
%                     i1, sprois{i2}, ...
%                     med_PWS, diffDir, med_PFS, p_rs(i1,i2));
        end
    end
end
fprintf(1, '\n');

return
            
