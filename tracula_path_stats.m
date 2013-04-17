function path_stats = tracula_path_stats(TRACULA_DIR, subjID, hemi, tractName, varargin)
%% Config: paths
% TRACULA_DIR = '/users/cais/STUT/tracula/';
check_dir(TRACULA_DIR);

%%
tract_dir = fullfile(TRACULA_DIR, subjID, 'dpath', [hemi, '.', tractName, '_*']);

d1 = dir(tract_dir);
if ~isempty(d1)
    tract_dir = fullfile(TRACULA_DIR, subjID, 'dpath', d1(1).name);
else
    error('Tract directory not found.');
end

stats_fn = fullfile(tract_dir, 'pathstats.overall.txt');
fprintf('stats_fn = %s\n', stats_fn);

if ~isfile(stats_fn)
    error('By-voxel stats file not found.');
end

strs = textread(stats_fn, '%s', 'delimiter', '\n');

for i1 = 1 : numel(strs)
    if isequal(strs{i1}(1), 'x')
        break;
    end
end

path_stats = struct;

for i1 = 1 : numel(strs);
    str = strs{i1};
    ss = splitstring(str);
    
    if isequal(deblank(ss{1}), 'Volume')
        path_stats.volume = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'Len_Min')
        path_stats.Len_Min = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'Len_Max')
        path_stats.Len_Max = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'Len_Avg')
        path_stats.Len_Avg = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'Len_Center')
        path_stats.Len_Center = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'AD_Avg')
        path_stats.AD_Avg = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'RD_Avg')
        path_stats.RD_Avg = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'FA_Avg')
        path_stats.FA_Avg = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'FA_Avg_Weight')
        path_stats.FA_Avg_Weight = str2num(ss{2});
    elseif isequal(deblank(ss{1}), 'FA_Avg_Center')
        path_stats.FA_Avg_Center = str2num(ss{2});
    end
end

return