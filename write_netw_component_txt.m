function [nEdges, nNodes] = ...
    write_netw_component_txt(adj, ncomp, rois, pvals, txtfn)
txtf = fopen(txtfn, 'wt');

assert(isequal(size(adj), size(pvals)));

idx = find(adj == ncomp);
% nEdges = length(idx);

N = numel(rois);

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
    
    fprintf(txtf, '%s - %s: sig=%f\n', ...
            rois{ix}, rois{iy}, -log10(t_pval));

    if isempty(fsic(nodeList, rois{ix}))
        nodeList{end + 1} = rois{ix};
    end
    if isempty(fsic(nodeList, rois{iy}))
        nodeList{end + 1} = rois{iy};
    end
    
	nEdges = nEdges + 1;
end

nNodes = length(nodeList);
fclose(txtf);
check_file(txtfn);

return