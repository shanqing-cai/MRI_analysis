function tract_subcort_conn(scSeed, mainMaskRatio, varargin)
%% Constants
tractSubcortConnDir = '/users/cais/STUT/analysis/tract_subcort_conn';
tractSegDir = '/users/cais/STUT/analysis/tractseg_aparc12';

grps = {'PWS', 'PFS'};

%% Config
P_THRESH_UNC = 0.05;

%%
cwd = pwd;
assert(isequal(cwd(end - 6 : end), 'scripts'));
check_dir(tractSubcortConnDir);

seedDir = fullfile(tractSubcortConnDir, scSeed);
check_dir(seedDir);

% -- Load the list of cortical ROIs -- %
rois = aparc12_cortical_rois();

nPerm = 0;
if ~isempty(fsic(varargin, '--perm'))
    nPerm = varargin{fsic(varargin, '--perm') + 1};
    assert(nPerm > 0);
    
    check_dir('perm_files', '-create');
    
    permMatFN = fullfile('perm_files', ...
        sprintf('%s_%s_%.2f_perm%d.mat', ...
                mfilename, scSeed, ...
                mainMaskRatio, nPerm));
end

%% -- Figure out the subject IDs and group labels -- %
check_dir(tractSegDir);
ds = dir(fullfile(tractSegDir, 'S*'));

sIDs.PWS = {};
sIDs.PFS = {};
SSI4 = [];
for i1 = 1 : numel(ds)
    if length(ds(i1).name) ~= 3
        continue;
    end
    
    t_sID = LUT_MRI_subjIDs(ds(i1).name, 'inv');
    t_grp = t_sID(1 : 3);
    
    sIDs.(t_grp){end + 1} = ds(i1).name;
    
    if isequal(t_grp, 'PWS')
        SSI4(end + 1) = get_qdec_measure(ds(i1).name, 'SSI');
    end
end


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

sp_rhos = nan(1, numel(rois));
sp_ts = nan(1, numel(rois));
sp_ps = nan(1, numel(rois));

medDiffs = nan(1, numel(rois));

hiliteROIs = {};
hiliteROIClrs = {};
drawBndClrs = {};

rp_sigs = nan(1 + nPerm, numel(rois));

for i0 = 1 : 1 + nPerm
    dat_a = [sccConn.PWS; sccConn.PFS];
    
    if i0 == 1 % No permutation
        
    else % Perform permutation
        idxperm = randperm(size(dat_a, 1));
        dat_a = dat_a(idxperm, :);
    end
    dat_PWS = dat_a(1 : size(sccConn.PWS, 1), :);
    dat_PFS = dat_a(size(sccConn.PWS, 1) + 1 : end, :);
        
    for i1 = 1 : numel(rois)
        roi = rois{i1};

        [rs_ps(i1), ~, stats] = ranksum(dat_PWS(:, i1), dat_PFS(:, i1));
        rs_zs(i1) = stats.zval;
        rs_rs(i1) = stats.ranksum;

        medDiffs(i1) = median(dat_PWS(:, i1)) - median(dat_PFS(:, i1));
        rs_sigs(i1) = sign(medDiffs(i1)) * -log10(rs_ps(i1));

        if i0 == 1
            if rs_ps(i1) < P_THRESH_UNC
                fprintf(1, '%s: ranksum z=%.3f; p=%.3f\n\t(Median: PWS=%f; PFS=%f)\n\t(IQR; PWS=%f; PFS=%f)\n', ...
                        roi, rs_zs(i1), rs_ps(i1), ...
                        median(dat_PWS(:, i1)), median(dat_PFS(:, i1)), ...
                        iqr(dat_PWS(:, i1)), iqr(dat_PFS(:, i1)));

                hiliteROIs{end + 1} = roi;
                hiliteROIClrs{end + 1} = [0, 1, 0];
                drawBndClrs{end + 1} = [0, 1, 0];
            end
        end

        % -- Correlation with SSI4 -- %
        [sp_rhos(i1), sp_ts(i1), sp_ps(i1)] = spear(SSI4(:), dat_PWS(:, i1));

        rp_sigs(i0, i1) = rs_sigs(i1);
    %     t_fillROIs{end + 1} = roi;
    %     t_fillClrs{end + 1} = 
    end
end

%%
if nPerm > 0
    save(permMatFN, 'rp_sigs');
    check_file(permMatFN);
end

nrp_sigs = rp_sigs(1, :);
rp_sigs = rp_sigs(2 : end, :);
max_rp_sigs = max(abs(rp_sigs), [], 2);

%%
% -- Print significant correlation results -- %
hashROIs = {};
hashSigns = [];
fprintf(1, '=== Results of SSI4 correlation ===\n');
for i1 = 1 : numel(rois)
    roi = rois{i1};
    if sp_ps(i1) < P_THRESH_UNC
        fprintf(1, '%s: spear: rho=%.3f; p=%.3f\n', ...
                roi, sp_rhos(i1), sp_ps(i1));
            
        hashROIs{end + 1} = roi;
        hashSigns(end + 1) = sign(sp_rhos(i1));
    end
end

% [TODO: Implement random permutation test]

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
                             'drawBndROIs', hiliteROIs, drawBndClrs, ...
                             'hashROIs', hashROIs, hashSigns);

figure;
imshow(imFASigMap.(hemi));

return