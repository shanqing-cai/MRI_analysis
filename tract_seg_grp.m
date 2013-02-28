function tract_seg_grp(scROI)
%% CONFIG: Subjects
sIDs.PWS = {'S01', 'S04', 'S06', 'S07', 'S08', 'S09', 'S10', 'S12', 'S15',  ...
            'S16', 'S17', 'S20', 'S21', 'S26', 'S28', 'S29', 'S33', 'S34', ...
            'S36', 'S37'};
            % Included S16 (But S16 should probably be kept in the group FA comparison) 
            % What's wrong with S21? Why was he excluded before, in
            % aparcSL_FA_analysis?
sIDs.PFS = {'S02', 'S03', 'S05', 'S11', 'S13', 'S14', 'S18', 'S19', 'S22', ...
            'S23', 'S25', 'S27', 'S30', 'S31', 'S32', 'S35', 'S39'};

%% CONFIG: Directories
tractSegBase = '/users/cais/STUT/analysis/tractseg_aparc12';

%% CONFIG: Misc.
pctErrTol = 0.1;

colors.PWS = [1, 0, 0];
colors.PFS = [0, 0, 1];

%%
grps = fields(sIDs);

lobeNames = {'Prefrontal', 'Premotor', 'Precentral', 'Postcentral', ...
             'PPC', 'Occipital', 'Temporal'};
         
nLobes = length(lobeNames);
lobeVolPct = struct;

for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    lobeVolPct.(grp) = nan(length(sIDs.(grp)), nLobes);
    
    for i2 = 1 : numel(sIDs.(grp))
        sID = sIDs.(grp){i2};
        
        sumTxtFN = fullfile(tractSegBase, sID, scROI, ...
                            'tract_regionParc.sum.txt');
                        
        if ~isfile(sumTxtFN)
            error('Cannot find summary file for subject %s: %s', sID, sumTxtFN);
        end
        
        txt = textread(sumTxtFN, '%s', 'delimiter', '\n');
        
        t_volPct = nan(1, nLobes);
        % --- AD HOC warning! --- %
        for j1 = 1 : nLobes            
            t_items = splitstring(txt{j1 + 2});
            t_volPct(j1) = str2double(t_items{end - 1});
        end
        % --- ~AD HOC warning! --- %
        
        % --- Sanity check --- %
        if (sum(t_volPct) < 100 - pctErrTol) || (sum(t_volPct) > 100 + pctErrTol)
            error('Percentages do not sum to 100 in subject %s', sID);
        end
        
        lobeVolPct.(grp)(i2, :) = t_volPct;
    end
end

%% Visualization
figure('Color', 'w');
for i1 = 1 : 2
    grp = grps{i1};
    semilogy(mean(lobeVolPct.(grp)), 'o-', 'color', colors.(grp));
    hold on;
end
grid on;
set(gca, 'XTick', [1 : 7], ...
    'XTickLabel', {'Prefrontal', 'Premotor', 'Precentral', 'Postcentral', 'PPC', 'Occipital', 'Temporal'});
xlabel('Predominant target of projection');
ylabel('Fraction of Putamen volume');
legend(grps);
    

return