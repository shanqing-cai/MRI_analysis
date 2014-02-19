function roi_names = get_aparc12_cortical_rois(varargin)
roi_names = {'adPMC', 'pMFg', 'pIFs', 'dIFo', 'vIFo', 'pFO', 'aINS', 'mdPMC', ...
        'pdPMC', 'midPMC', 'vPMC', 'aCO', 'pINS', 'dMC', 'vMC', ...
        'pCO', 'dSC', 'vSC', 'aSMg', 'PO', ...
        'TP', 'PP', 'Hg', 'PT', 'pSTg', 'pdSTs', 'adSTs', 'pvSTs', 'avSTs', ...
        'aMTg', 'pMTg', ...
        'SMA', 'preSMA', 'aCGg', ...
        'FP', 'aMFg', 'aIFs', 'SFg', 'aFO', 'FOC', 'SPL', 'pSMg', ...
        'Ag', 'OC', 'MTOg', 'ITOg', 'pITg', 'PCN', 'Lg', 'pPH', ...
        'aPH', 'midCGg', 'pCGg', 'FMC', ...
        'midMC', 'aSTg'}; % Added (2013-05-07)
    
% speech_roi_names_old = {'Hg', 'PO', 'PP', 'PT', 'SMA', ...
%                     'aCGg', 'aCO', 'aFO', 'aINS', 'aSMg', ...
%                     'dIFo', 'dMC', 'dSC', 'midMC', 'midPMC', ...
%                     'pCO', 'pFO', 'pIFs', 'pINS', 'aSTg', 'pSTg', ...
%                     'pdPMC', 'pdSTs', 'preSMA', 'vIFo', 'vMC', ...
%                     'vPMC', 'vSC'};
% --- Ordered by Lobe and AP position --- %
speech_roi_names = {'aFO', 'pFO', 'vIFo', 'dIFo', 'pIFs', ...
                    'vPMC', 'midPMC', 'preSMA', 'SMA', 'pdPMC', ...
                    'aINS', 'pINS', ...
                    'aCO', 'vMC', 'midMC', 'dMC', ...
                    'pCO', 'vSC', 'dSC', 'PO', 'aSMg', ...
                    'aSTg', 'PP', 'pdSTs', 'Hg', 'pSTg', 'PT', 'aCGg'};
                
% -- DEBUG -- %
% sns = sort(speech_roi_names);
% sns_old = sort(speech_roi_names_old);
% for i1 = 1 : numel(sns)
%     if ~isequal(sns{i1}, sns_old{i1})
%         pause(0);
%     end
% end
% -- ~DEBUG -- %
% In addition: 'aSTg'

stutFL_roiSet_fn_wc = '/users/cais/STUT/scripts/activROIs_uc_thr%.3f.mat';


if nargin >= 1 && isequal(varargin{1}, 'speech')
    roi_names = speech_roi_names;
elseif nargin >= 1 && ...
       (isequal(varargin{1}(1 : 13), 'speech_PFS_lh') || ...
        isequal(varargin{1}(1 : 13), 'speech_PFS_rh') || ...
        isequal(varargin{1}(1 : 13), 'speech_PWS_lh') || ...
        isequal(varargin{1}(1 : 13), 'speech_PWS_rh') || ...
        isequal(varargin{1}(1 : 12), 'speech_2g_lh') || ...
        isequal(varargin{1}(1 : 12), 'speech_2g_rh'))
    idxus = strfind(varargin{1}, '_');
    t_grp = varargin{1}(idxus(1) + 1 : idxus(2) - 1);
    t_hemi = varargin{1}(idxus(2) + 1 : idxus(3) - 1);
    t_thr = str2double(varargin{1}(idxus(3) + 1 : end));
    
    stutFL_roiSet_fn = sprintf(stutFL_roiSet_fn_wc, t_thr);
    
    load(stutFL_roiSet_fn); % gives roiSet;
    if ~isequal(t_grp, '2g')
        roi_names = sort(roiSet.(t_grp).(t_hemi));
    else
        roi_names_PFS = sort(roiSet.PFS.(t_hemi));
        roi_names_PWS = sort(roiSet.PWS.(t_hemi));
        roi_names = union(roi_names_PFS, roi_names_PWS);      
    end
end

if ~isempty(fsic(varargin, 'lh'))
    hemi = 'lh';
elseif ~isempty(fsic(varargin, 'rh'))
    hemi = 'rh';
elseif ~isempty(fsic(varargin, 'xh'))
    hemi = 'xh';
else
    hemi = '';
end

if ~isempty(hemi) && ~isequal(hemi, 'xh')
    for i1 = 1 : numel(roi_names)
        roi_names{i1} = [hemi, '_', roi_names{i1}];
    end
else
    hemis = {'lh', 'rh'};
    nr = length(roi_names);
    roi_names = [roi_names, roi_names];
    for i0 = 1 : numel(hemis)
        t_hemi = hemis{i0};
        for i1 = 1 : nr
            roi_names{(i0 - 1) * nr + i1} = [t_hemi, '_', roi_names{(i0 - 1) * nr + i1}];
        end
    end
end


return