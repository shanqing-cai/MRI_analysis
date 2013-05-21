function bMatch = roi_name_match(rn1, rn2)
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
it1 = splitstring(rn1, '_');
hemi1 = it1{1};
name1 = it1{2};

it2 = splitstring(rn2, '_');
hemi2 = it2{1};
name2 = it2{2};

hemiMatch = isequal(hemi1, hemi2);

if ~isempty(fsic(ali(:, 1), name1))
    idx = fsic(ali(:, 1), name1);
    
    nameMatch = isequal(name1, name2) || isequal(ali{idx, 2}, name2);
    
elseif ~isempty(fsic(ali(:, 2), name1))
    idx = fsic(ali(:, 2), name1);
    
    nameMatch = isequal(name1, name2) || isequal(ali{idx, 1}, name2);
    
else
    nameMatch = isequal(name1, name2)
end

bMatch = hemiMatch & nameMatch;

return