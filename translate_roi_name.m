function rn1 = translate_roi_name(rn0)
%%
ali = {'H',     'Hg'; ...
       'aCG',   'aCGg'; ...
       'midCG', 'midCGg'; ...
       'pCG',   'pCGg'; ...
       'LG',    'Lg'; ...
       'AG',    'Ag'; ...
       'ITO',   'ITOg'; ...
       'MTO',   'MTOg'};

%%
it = splitstring(rn0, '_');

if length(it) == 2
    hemi = it{1};
    rn = it{2};
elseif length(it) == 1
    rn = it{1};
else
    error('Unrecognized format in ROI name: %s', rn0);
end
    
if ~isempty(fsic(ali(:, 1), rn))
    rn1 = ali{fsic(ali(:, 1), rn), 2};
else
    rn1 = rn;
end

if length(it) == 2
    rn1 = [hemi, '_', rn1];
else
    rn1 = rn1;
end
return