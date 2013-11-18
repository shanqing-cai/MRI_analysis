function [nEdges, nNodes] = ...
    write_netw_component_txt(adj, ncomp, rois, pvals, txtfn, varargin)
%% Process optional input arguments 
bXH = ~isempty(fsic(varargin, '--xh'));

%%
txtf = fopen(txtfn, 'wt');

if bXH && (size(pvals, 1) == size(adj, 1) / 2) ...
       && (size(pvals, 2) == size(adj, 2) / 2)
    pvals0 = pvals;
    N0 = size(pvals, 1);
    pvals = zeros(N0 * 2, N0 * 2);
    pvals(1 : N0, N0 + 1 : 2 * N0) = pvals0;
    pvals = pvals + pvals';
end
assert(isequal(size(adj), size(pvals)));

idx = find(adj == ncomp);
% nEdges = length(idx);

% N = numel(rois);
N = size(adj, 1);

nEdges = 0;
nNodes = 0;
nodeList = {};

for i1 = 1 : length(idx)    
    [ix, iy] = ind2sub([N, N], idx(i1));        
    
    if ix >= iy
        continue;
    end
    
    if isnan(pvals(ix, iy))
        t_pval = pvals(iy, ix);
        if isnan(t_pval)
            error('Unable to find non-NaN p-value');
        end
    else
        t_pval = pvals(ix, iy);
    end
    
    if ~bXH    
        fprintf(txtf, '%s - %s: sig=%f\n', ...
                rois{ix}, rois{iy}, -log10(t_pval));
    else
        fprintf(txtf, '%s - %s: sig=%f\n', ...
                sprintf('lh_%s', rois{ix}), ...
                sprintf('rh_%s', rois{iy - N0}), ...
                -log10(t_pval));
    end

    if isempty(fsic(nodeList, rois{ix}))
        nodeList{end + 1} = rois{ix};
    end
    
    if ~bXH
        iy1 = iy;
    else
        iy1 = iy - N0;
    end
    
    if isempty(fsic(nodeList, rois{iy1}))
        nodeList{end + 1} = rois{iy1};
    end
    
	nEdges = nEdges + 1;
end

nNodes = length(nodeList);
fclose(txtf);
check_file(txtfn);

return