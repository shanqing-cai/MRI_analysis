function show_tract_regionSeg_3d(inMatFN, alpha)
%% CONFIG
common_graph_path = 'E:\speechres\commonmcode\graph';

% alpha = 1;

%%
paths = path;
if isempty(strfind(paths, common_graph_path))
    addpath(common_graph_path);
end

%%
load(inMatFN); % gives vol, cdata

vol = flipdim(vol, 1);

alphaMat = alpha * ones(size(vol));
alphaMat(vol == 0) = 0;

vol3d('cdata', vol, 'alpha', alphaMat);
axis equal;
grid on;

colorbar;
return