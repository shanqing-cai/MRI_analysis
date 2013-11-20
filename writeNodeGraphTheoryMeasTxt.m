function writeNodeGraphTheoryMeasTxt(txtFN, rois, ...
                                     a_str, p_str, ...
                                     a_bc, p_bc, ...
                                     a_cc, p_cc)
%%
% Input arguments: txtFN: text file namt
%                  rois: list of ROIs
%                  a_str, p_str: node strenghts and BGC p-values,
%                  a_bc, p_bc: node betweenness centrality (BC) and BGC
%                  p-values
%                  a_cc, p_cc: node clustering coefficient (CC) and BGC
%                  p-values

%% Configure
measures = {'str', 'Node strengths',                '%.1f'; ...
            'bc',  'Node betweenness centrality',   '%.4f'; ...
            'cc',  'Node clustering coefficient',   '%.2f'};

%% 
rois1 = cell(1, numel(rois));
for i1 = 1 : numel(rois)
    rois1{i1} = strrep(strrep(rois{i1}, 'lh_', 'L '), 'rh_', 'R ');
end

%%
f = fopen(txtFN, 'wt');


grps = fields(a_str);
grps = sort(grps);
        
for i1 = 1 : size(measures, 1)
    eval(sprintf('a_x = a_%s;', measures{i1, 1}));
    eval(sprintf('ps = p_%s;', measures{i1, 1}));
    fmt = measures{i1, 3};
    
    assert(length(rois1) == length(ps));
    
    mn_x = struct;
    sd_x = struct;
    
    for i2 = 1 : numel(grps)
        grp = grps{i2};
        
        assert(length(rois1) == size(a_x.(grp), 1));
        
        mn_x.(grp) = nanmean(a_x.(grp), 2);
        sd_x.(grp) = nan(size(mn_x.(grp)));
        
        for i3 = 1 : size(mn_x.(grp), 1)
            sd_x.(grp)(i3) = nanstd(a_x.(grp)(i3, :));
        end
    end
    
    fprintf(f, '=== %s ===\n', measures{i1, 2});
    fprintf(f, 'Node\t\t');
    for i2 = 1 : numel(grps)
        grp = grps{i2};
        fprintf(f, '%s\t\t', grp);
    end
    fprintf(f, 'Significance');
    fprintf(f, '\n');
    
    for i2 = 1 : numel(rois1)        
        fprintf(f, '%s\t', rois1{i2});
        
        for i3 = 1 : numel(grps)
            grp = grps{i3};
            
            fprintf(f, [fmt, '+/-', fmt, '\t'], mn_x.(grp)(i2), sd_x.(grp)(i2));
        end
        
        %--- Plot significance ---%
        if ps(i2) < 0.005
            fprintf(f, '***');
        elseif ps(i2) < 0.01
            fprintf(f, '**');
        elseif ps(i2) < 0.05
            fprintf(f, '*');
        end
        
        fprintf(f, '\n');
    end
    
    
    
    fprintf(f, '\n');
end

fclose(f);
%%
fprintf(1, 'Written to file %s\n', txtFN);

return