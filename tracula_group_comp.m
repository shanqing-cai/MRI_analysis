function tracula_group_comp(avgFAMode)
% Input arguments:
%    argFAMode: weight or center
% 
%% Config: subject IDs   
SIDS.PWS = {'S01', 'S04', 'S06', 'S07', 'S08', 'S09', 'S10', 'S12', ...
            'S15', 'S16', 'S17', 'S20', 'S21', 'S26', 'S28', 'S29', ...
            'S33', 'S34', 'S36', 'S37'};
% SIDS.PWS = {'S01', 'S04', 'S06', 'S07', 'S08', 'S09', 'S10', 'S12', ...
%             'S15', 'S16', 'S17', 'S20', 'S21', 'S26', 'S28', 'S29'};
SIDS.PFS = {'S02', 'S03', 'S05', 'S11', 'S13', 'S14', 'S18', 'S19', ...
            'S22', 'S23', 'S25', 'S27', 'S30', 'S31', 'S32', ...
            'S35', 'S39'}; % Left out S24; S38

TRACULA_RES_DIR = '/users/cais/STUT/analysis/dti2/tracula';

b_500 = 0;
b_ns1e4 = 0;

TRACTS = {'lh_slft', 'lh_slfp', 'lh_cst', 'lh_ilf', 'lh_atr', 'lh_unc', ...
          'rh_slft', 'rh_slfp', 'rh_cst', 'rh_ilf', 'rh_atr', 'rh_unc'};

%% Config: data processing options
INTERP_N = 200;
VOXEL_VOL = 2 * 2 * 2; % mm^3

%% Version 5.0.0 or later?
check_dir(TRACULA_RES_DIR);

grps = fields(SIDS);
if b_500 == 1
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        for i2 = 1 : numel(SIDS.(grp))
            SIDS.(grp){i2} = [SIDS.(grp){i2}, '_500'];
        end
    end
end
if b_ns1e4 == 1
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        for i2 = 1 : numel(SIDS.(grp))
            SIDS.(grp){i2} = [SIDS.(grp){i2}, '_ns1e4_1'];
        end
    end
end

%%
protoArray = struct;
protoArrayN = struct;


for i1 = 1 : numel(grps)
    grp = grps{i1};
    protoArray.(grp) = nan(1, numel(SIDS.(grp)));
    protoArrayN.(grp) = nan(numel(SIDS.(grp)), INTERP_N);
end

%%
avgTractFA = struct;
tractVol = struct;
tractFA_bv = struct;    % By-voxel FA (interpolated)
tractLenAvg = struct; 

for i0 = 1 : numel(TRACTS)
    t_tract = TRACTS{i0};
    idx_hf = strfind(t_tract, '_');
    idx_hf = idx_hf(1);
    hemi = t_tract(1 : idx_hf(1) - 1);
    roi = t_tract(idx_hf(1) + 1 : end);
    
    avgTractFA.(t_tract) = protoArray;
    tractVol.(t_tract) = protoArray;
    tractFA_bv.(t_tract) = protoArrayN;
    
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        for i2 = 1 : numel(SIDS.(grp))
            path_stats = ...
                tracula_path_stats(TRACULA_RES_DIR, SIDS.(grp){i2}, hemi, roi);
            [tr_x, tr_y, tr_z, tr_FA, norm_dist, norm_tr_FA] = ...
                tracula_path_bv(TRACULA_RES_DIR, SIDS.(grp){i2}, hemi, roi);
            % tr_FA = interp1(1 : numel(tr_FA), tr_FA, linspace(1, numel(tr_FA), INTERP_N));            
            % tractFA_bv.(t_tract).(grp)(i2, :) = tr_FA;
            tractFA_bv.(t_tract).(grp)(i2, :) = norm_tr_FA;
            
            if isequal(avgFAMode, 'weight')
                avgTractFA.(t_tract).(grp)(i2) = path_stats.FA_Avg_Weight;
            elseif isequal(avgFAMode, 'center')
                avgTractFA.(t_tract).(grp)(i2) = path_stats.FA_Avg_Center;
            else
                error('Unrecognized avgFAMode: %s', avgFAMode);
            end
            
            tractVol.(t_tract).(grp)(i2) = path_stats.volume * VOXEL_VOL;
        end
    end
end

% Calculate the ICV (intra-cranial volume) and STV (supra-tentorial volume)eicv = struct;
stv = struct;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    icv.(grp) = nan(1, numel(SIDS.(grp)));
    stv.(grp) = nan(1, numel(SIDS.(grp)));
    
    for i2 = 1 : numel(SIDS.(grp))
        [t_icv, t_stv] = get_ICV_STV(SIDS.(grp){i2});
        icv.(grp)(i2) = t_icv;
        stv.(grp)(i2) = t_stv;
    end
end

% Calculate normalized tract volumes (normalized by ICV)
tractVol_normICV = struct;
trcts = fields(tractVol);
for i0 = 1 : numel(trcts)
    trct = trcts{i0};
    tractVol_normICV.(trct) = struct;
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        tractVol_normICV.(trct).(grp) = tractVol.(trct).(grp) ./ icv.(grp);
    end
end

%% Visualization: tract avg. FAs
figure('Position', [50, 50, 1000, 600]);
subplot(2, 3, 1); meta_comp2(avgTractFA.lh_cst, 'L. CST avg. FA', 'noFig');
subplot(2, 3, 2); meta_comp2(avgTractFA.lh_slft, 'L. SLFT avg. FA', 'noFig');
subplot(2, 3, 3); meta_comp2(avgTractFA.lh_slfp, 'L. SLFP avg. FA', 'noFig');
subplot(2, 3, 4); meta_comp2(avgTractFA.lh_ilf, 'L. ILF avg. FA', 'noFig');
subplot(2, 3, 5); meta_comp2(avgTractFA.lh_atr, 'L. ATR avg. FA', 'noFig');
subplot(2, 3, 6); meta_comp2(avgTractFA.lh_unc, 'L. UNC avg. FA', 'noFig');

figure('Position', [50, 50, 1000, 600]);
subplot(2, 3, 1); meta_comp2(avgTractFA.rh_cst, 'R. CST avg. FA', 'noFig');
subplot(2, 3, 2); meta_comp2(avgTractFA.rh_slft, 'R. SLFT avg. FA', 'noFig');
subplot(2, 3, 3); meta_comp2(avgTractFA.rh_slfp, 'R. SLFP avg. FA', 'noFig');
subplot(2, 3, 4); meta_comp2(avgTractFA.rh_ilf, 'R. ILF avg. FA', 'noFig');
subplot(2, 3, 5); meta_comp2(avgTractFA.rh_atr, 'R. ATR avg. FA', 'noFig');
subplot(2, 3, 6); meta_comp2(avgTractFA.rh_unc, 'R. UNC avg. FA', 'noFig');

% -- Special figure: avg. FA of the R SLFT -- 
figure('Position', [100, 100, 300, 300]);
meta_comp2(avgTractFA.rh_slft, 'Avg. FA in R arcuate fasciculus', ...
           'showMeanSE', 'PFS_first', 'noFig', 'FontSize', 12);


%% Visualization: tract volume
figure('Position', [50, 50, 800, 600]);
subplot(2, 2, 1); meta_comp2(tractVol.lh_cst, 'L. CST volume', 'noFig');
subplot(2, 2, 2); meta_comp2(tractVol.lh_slft, 'L. SLFT volume', 'noFig');
subplot(2, 2, 3); meta_comp2(tractVol.lh_slfp, 'L. SLFP volume', 'noFig');
subplot(2, 2, 4); meta_comp2(tractVol.lh_ilf, 'L. ILF volume', 'noFig');

figure('Position', [50, 50, 800, 600]);
subplot(2, 2, 1); meta_comp2(tractVol.rh_cst, 'R. CST volume', 'noFig');
subplot(2, 2, 2); meta_comp2(tractVol.rh_slft, 'R. SLFT volume', 'noFig');
subplot(2, 2, 3); meta_comp2(tractVol.rh_slfp, 'R. SLFP volume', 'noFig');
subplot(2, 2, 4); meta_comp2(tractVol.rh_ilf, 'R. ILF volume', 'noFig');

% -- Special figure: volume of the R SLFT -- 
figure('Position', [100, 100, 300, 300]);
meta_comp2(tractVol.rh_slft, 'Volume of R arcuate fasciculus', ...
           'showMeanSE', 'PFS_first', 'noFig', 'FontSize', 12);

%% Visualization: tract volume normalized by ICV
figure('Position', [60, 60, 750, 500]);
subplot(2, 2, 1); meta_comp2(tractVol_normICV.lh_slft, 'L. SLFT volume (norm. by ICV)', 'noFig');
subplot(2, 2, 2); meta_comp2(tractVol_normICV.rh_slft, 'R. SLFT volume (norm. by ICV)', 'noFig');
subplot(2, 2, 3); meta_comp2(tractVol_normICV.lh_ilf, 'L. ILF volume (norm. by ICV)', 'noFig');
subplot(2, 2, 4); meta_comp2(tractVol_normICV.rh_ilf, 'R. ILF volume (norm. by ICV)', 'noFig');

%% Correlation: SSI-4 and avg FA
SSI4.PWS = ds_SSI4(SIDS.PWS);

figure('Position', [70, 70, 600, 400]);
subplot(2, 2, 1); meta_corr(SSI4.PWS, [], avgTractFA.rh_slft.PWS, [], 'SSI4', 'R. SLFT avg. tract FA', 'noFig');
subplot(2, 2, 2); meta_corr(SSI4.PWS, [], tractVol.rh_slft.PWS, [], 'SSI4', 'R. SLFT volume', 'noFig');
subplot(2, 2, 3); meta_corr(SSI4.PWS, [], avgTractFA.lh_ilf.PWS, [], 'SSI4', 'L. ILF avg. tract FA', 'noFig');
%% Correaltion: with EH_comp
% EH_comp.PWS = ds_EH_comp(SIDS.PWS);
% EH_comp.PFS = ds_EH_comp(SIDS.PFS);
EH_comp.PWS = get_qdec_measure(SIDS.PWS, 'EH_comp_300');
EH_comp.PFS = get_qdec_measure(SIDS.PFS, 'EH_comp_300');

meta_corr(EH_comp.PWS, EH_comp.PFS, avgTractFA.rh_slft.PWS, avgTractFA.rh_slft.PFS, ...
          'EH_comp', 'avg. tract FA');
meta_corr(EH_comp.PWS, [], avgTractFA.rh_slft.PWS, [], ...
          'EH_comp', 'avg. tract FA');
      
%% Correlation: with T comp (U1 timing)
% T_comp.PWS = ds_T_comp(SIDS.PWS, 'IUIYInt_decelAccel');
% T_comp.PFS = ds_T_comp(SIDS.PFS, 'IUIYInt_decelAccel');
% T_comp.PWS = ds_T_comp(SIDS.PWS, 'IUInt_decelAccel');
% T_comp.PFS = ds_T_comp(SIDS.PFS, 'IUInt_decelAccel');
T_comp.PWS = get_qdec_measure(SIDS.PWS, 'tempResp');
T_comp.PFS = get_qdec_measure(SIDS.PFS, 'tempResp');

meta_corr(T_comp.PWS, T_comp.PFS, avgTractFA.rh_slft.PWS, avgTractFA.rh_slft.PFS, ...
          'Temporal perturbation response', 'R. SLFT avg. tract FA'); 
meta_corr(T_comp.PWS, T_comp.PFS, tractVol.rh_slft.PWS, tractVol.rh_slft.PFS, ...
          'Temporal perturbation response', 'R. SLFT tract volume');
      
%% Correlation: with T comp (U1Y1 avg)
T_comp_UY.PWS = get_qdec_measure(SIDS.PWS, 'tempRspUY');
T_comp_UY.PFS = get_qdec_measure(SIDS.PFS, 'tempRspUY');

meta_corr(T_comp_UY.PWS, T_comp_UY.PFS, avgTractFA.rh_slft.PWS, avgTractFA.rh_slft.PFS, ...
          'Temporal perturbation response (U1Y1 avg)', 'R. SLFT avg. tract FA'); 
meta_corr(T_comp_UY.PWS, T_comp_UY.PFS, tractVol.rh_slft.PWS, tractVol.rh_slft.PFS, ...
          'Temporal perturbation response (U1Y1 avg)', 'R. SLFT tract volume');

%% Visualization: by-voxel FA
meta_bv_fa2(tractFA_bv.lh_slft, 'L. SLFT by-voxel FA', 'corr', T_comp);
meta_bv_fa2(tractFA_bv.rh_slft, 'R. SLFT by-voxel FA', 'corr', T_comp);

meta_bv_fa2(tractFA_bv.lh_cst, 'L. CST by-voxel FA', 'corr', T_comp);
meta_bv_fa2(tractFA_bv.rh_cst, 'R. CST by-voxel FA', 'corr', T_comp);

return

%%
function meta_corr(x_PWS, x_PFS, y_PWS, y_PFS, measName1, measName2, varargin)
if isempty(fsic(varargin, 'noFig'))
    figure('Color', 'w');
end
plot(x_PWS, y_PWS, 'ro');
hold on;
plot(x_PFS, y_PFS, 'ko');

x = [x_PWS, x_PFS];
y = [y_PWS, y_PFS];
[k, r2, p] = lincorr(x, y);
if ~isempty(x_PWS)
    [k_PWS, r2_PWS, p_PWS] = lincorr(x_PWS, y_PWS);
else
    k_PWS = NaN;
    r2_PWS = NaN;
    p_PWS = NaN;
end
if ~isempty(x_PFS)
    [k_PFS, r2_PFS, p_PFS] = lincorr(x_PFS, y_PFS);
else
    k_PFS = NaN;
    r2_PFS = NaN;
    p_PFS = NaN;
end

xlabel(measName1); 
ylabel(measName2);

xs = get(gca, 'XLim');
ys = get(gca, 'YLim');
plot(xs, xs * k(2) + k(1), '--');

text(xs(1) + 0.1 * range(xs), ys(2) - 0.1 * range(ys), ...
     sprintf('Two-group Lin-corr: p = %.4f', p));
return

