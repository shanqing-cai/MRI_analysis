function stutFL_gen_roi_timecourses_aparc12(sID, runNum)
%% CONFIG
% roiDataDir = '/users/cais/STUT/analysis/maclearn';
roiDataDir = '/users/cais/STUT/analysis/aparc12_FL';
funclocDir = '/users/cais/STUT/funcloc';
sched_fn = fullfile('/users/cais/STUT/DATA/', sID, 'funcloc_sched.txt');

pThresh_uc = 0.05;

%%
check_file(sched_fn);

sDir = fullfile(roiDataDir, sID);
% roiDataTxt = fullfile(sDir, 'r843000-11.txt');
roiDataWC = dir(fullfile(sDir, 'r*-*.txt'));
if runNum < 1 || runNum > numel(roiDataWC)
    error('Erroneous runNum');
end
roiDataTxt = fullfile(sDir, roiDataWC(runNum).name);
check_file(roiDataTxt);
fprintf('roiDataTxt = %s\n', roiDataTxt);

if ~isfile(sched_fn)
    error('Stim schedule file not found: %s', sched_fn);
end

[roiNames, roiMat] = read_roi_data(roiDataTxt);

n = size(roiMat, 1);

roi_timecourses = struct;
roi_timecourses.nROIs = numel(roiNames);
roi_timecourses.roiNames = roiNames;


%% Read the funcloc schedule (to be later compared with the temporal EV)
txt = textread(sched_fn, '%s', 'delimiter', '\n');
if isempty(strfind(txt{3 * (runNum - 1) + 1}, sprintf('Run %d: speech: ', runNum)))
    error('Unexpected format of sched file');
end
txt = strrep(txt{3 * (runNum - 1) + 1}, sprintf('Run %d: speech: ', runNum), '');
ttst = splitstring(txt); % timing of the speech trials
tst = nan(size(ttst, 1), 1);
for i1 = 1 : length(ttst)
    tst(i1) = str2double(ttst{i1});
end

stidc = zeros(1, n);  % speech-trial indicator
for i1 = 1 : numel(tst)
    stidc((tst(i1) - 2.5) / 5 + 1) = 1;
end

% stidc = stidc(1 : end - 1); % Empirically determined!!
stidc = stidc(2 : end);
nTrials_S = numel(find(stidc == 1));
nTrials_BL = numel(find(stidc == 0));

roi_timecourses.nFrames = numel(stidc);
roi_timecourses.dcv_data = nan(roi_timecourses.nROIs, roi_timecourses.nFrames);
roi_timecourses.dcv_data_S = nan(roi_timecourses.nROIs , nTrials_S);
roi_timecourses.dcv_data_BL = nan(roi_timecourses.nROIs , nTrials_BL);

% figplot(sig_avg); hold on; plot(stidc * 3, 'r');

% figure;
% plot(stidc(1 : end - 2), sig_l_vMC(3 : end), 'o-');

% figure;
% plot(stidc(1 : end - 1), sig_l_vMC(2 : end), 'o-');

%% Get words info
funclocWC = dir(fullfile(funclocDir, [sID, '*_*.mat']));
funclocFN = fullfile(funclocDir, funclocWC(runNum).name);
fprintf('funlocFN = %s\n', funclocFN);
load(funclocFN);    % gives data: 1 - 6 : {'Topic', 'Teacup','Boutique','######','######','######'};
widx = zeros(size(stidc));
widx(data.stimtype == 1) = 1;
widx(data.stimtype == 2) = 2;
widx(data.stimtype == 3) = 3;

%% Deconvolution
% 1. Construct the hrf matrix
% hrf= spm_hrf(1, [6, 16, 1, 1, 6, 0, 10])';
hrf_mat = nan(n - 1, n - 1);
TR = 5;

hrf_row = spm_hrf(TR, [6, 16, 1, 1, 6, - TR * 1.00, TR * (n - 1)])';
hrf_row(isnan(hrf_row)) = 0;

hrf_row = hrf_row(1 : end - 1);

for i1 = 1 : n - 1
    hrf_mat(i1, :) = [zeros(1, i1 - 1), hrf_row(1 : n - i1)];
end

% fit_l_vMC = stidc * hrf_mat;
% [k, r2, p] = lincorr(fit_l_vMC, sig_l_vMC);

% dcs_l_vMC = sig_l_vMC' * inv(hrf_mat);
nSigROIs = 0;
for i1 = 1 : roi_timecourses.nROIs
    t_sig = roiMat(:, i1);
    t_sig = t_sig(2 : end);
    t_sig = detrend(t_sig);
    dcs_sig = inv(hrf_mat) * t_sig;
    roi_timecourses.dcv_data(i1, :) = dcs_sig';
    roi_timecourses.dcv_data_S(i1, :) = dcs_sig(stidc == 1)';
    roi_timecourses.dcv_data_BL(i1, :) = dcs_sig(stidc == 0)';
    roi_timecourses.dcv_data_topic(i1, :) = dcs_sig(widx == 1)';
    roi_timecourses.dcv_data_teacup(i1, :) = dcs_sig(widx == 2)';
    roi_timecourses.dcv_data_boutique(i1, :) = dcs_sig(widx == 3)';
    
    [h, p] = ttest2(dcs_sig(stidc == 1), dcs_sig(stidc == 0));
    if p < pThresh_uc && mean(dcs_sig(stidc == 1)) > mean(dcs_sig(stidc == 0))
        fprintf('ROI %s: p = %f\n', roiNames{i1}, p);        
        nSigROIs = nSigROIs + 1;
    end
end
fprintf('%d of %d ROIs significant at p < %.3f\n', ...
        nSigROIs, roi_timecourses.nROIs, pThresh_uc)

%% Save to file
[fpath, fname] = fileparts(roiDataTxt);
idx_hf = strfind(fname, '-');
idx_dot = strfind(fname, '.');
realRunNum = str2num(fname(idx_hf + 1 : end));

saveFN = fullfile(sDir, ['roi_timecourses.', num2str(realRunNum), '.mat']);
save(saveFN, 'roi_timecourses');
fprintf('Data saved to %s\n', saveFN);

return