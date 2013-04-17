function [path_x, path_y, path_z, path_FA, norm_dist, norm_path_FA] = ...
            tracula_path_bv(TRACULA_DIR, subjID, hemi, tractName, varargin)
%% Config: paths
% TRACULA_DIR = '/users/cais/STUT/tracula/';
N_INTERP = 200;

%%
tract_dir = fullfile(TRACULA_DIR, subjID, 'dpath', [hemi, '.', tractName, '_*']);

d1 = dir(tract_dir);
if ~isempty(d1)
    tract_dir = fullfile(TRACULA_DIR, subjID, 'dpath', d1(1).name);
else
    error('Tract directory not found.');
end

bv_fn = fullfile(tract_dir, 'pathstats.byvoxel.txt');
fprintf('bv_fn = %s\n', bv_fn);

if ~isfile(bv_fn)
    error('By-voxel stats file not found.');
end

strs = textread(bv_fn, '%s', 'delimiter', '\n');

for i1 = 1 : numel(strs)
    if isequal(strs{i1}(1), 'x')
        break;
    end
end

strs = strs(i1 + 1 : end - 1);

path_x = [];
path_y = [];
path_z = [];
path_FA = [];
norm_dist = NaN;
norm_path_FA = NaN;

for i1 = 1 : numel(strs)
    str = strs{i1};
    ss = splitstring(str);
    path_x(end + 1) = str2num(ss{1});
    path_y(end + 1) = str2num(ss{2});
    path_z(end + 1) = str2num(ss{3});
    path_FA(end + 1) = str2num(ss{end});
end

%% Generate the FA along the normalized tract
d_x = diff(path_x); 
d_y = diff(path_y);
d_z = diff(path_z);

d_l = [0, cumsum(sqrt(d_x .^ 2 + d_y .^ 2 + d_z .^ 2))];
d_l = d_l / d_l(end);
norm_dist = linspace(0, 1, N_INTERP);  

if length(d_l) ~= length(path_FA)
    return
    error('Length mismatch.');
end

norm_path_FA = interp1(d_l, path_FA, norm_dist);

%%
bPlot = 0;
if ~isempty(fsic(varargin, 'plot'))
    bPlot = 1;
end

if bPlot
    figure;
    plot3(path_x, path_y, path_z);
    grid on;
    
    figure;
    plot(path_FA);
    ylabel('Fractional anisotropy');
end

return