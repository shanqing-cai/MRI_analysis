function [nSigs, sigConnections] = show_2d_mat(sig_rs, sprois, hemi, figName, visParams)
figSize = visParams(1);
verticalPadding = visParams(2); 
horizontalPadding = visParams(3);
cellShift = visParams(4);

%%
max_abs = max(abs(sig_rs(:)));
sig_rs = [sig_rs; zeros(1, size(sig_rs, 1))];
sig_rs(end, 1) = max_abs;
sig_rs(end, 2) = -max_abs;

figure('Position', [100, 300, figSize, figSize], 'Color', 'w', ...
       'Name', figName);

imagesc(sig_rs);

hold on;

% Create green-red colormap
cm = colormap;
ncm = size(cm, 1);
cm_upper = linspace(0, 1, ncm / 2 + 1);
% cm_upper = [cm_upper', ones(length(cm_upper), 2)];
cm_upper = [cm_upper', cm_upper', ones(length(cm_upper), 1)];
cm_lower = linspace(1, 0, ncm / 2 + 1);
% cm_lower = [ones(length(cm_lower), 1), cm_lower', ones(length(cm_lower), 1)];
cm_lower = [ones(length(cm_lower), 1), cm_lower', cm_lower'];
cm = [cm_upper(1 : end - 1, :); cm_lower(2 : end, :)];
colormap(cm);

sigConnections = {};
nSigs = [0, 0]; % [Pos, Neg]
for k1 = 1 : numel(sprois)
    for k2 = 1 : numel(sprois)
        if abs(sig_rs(k2, k1)) > abs(log10(0.005))
            text(k1 - cellShift, k2 - cellShift, '*', 'Color', 'w');
        elseif abs(sig_rs(k2, k1)) > abs(log10(0.01))
            text(k1 - cellShift, k2 - cellShift, 'X', 'Color', 'w');
        elseif abs(sig_rs(k2, k1)) > abs(log10(0.05))
            text(k1 - cellShift, k2 - cellShift, 'O', 'Color', 'w');        
        end
        
        if abs(sig_rs(k2, k1)) > abs(log10(0.05))
            if sig_rs(k2, k1) > 0
                nSigs(1) = nSigs(1) + 1;
            else
                nSigs(2) = nSigs(2) + 1;
            end
        end
        
        if abs(sig_rs(k2, k1)) > abs(log10(0.05))
            sigConnections{end + 1} = {sprois{k1}, sprois{k2}};
        end
    end
end

N = sum(nSigs);
chi2tab = [nSigs; [floor(N / 2), ceil(N / 2)]];
p_chi2 = chi2test(chi2tab);

set(gca, 'XLim', [0.5, length(sprois) + 0.5], ...
         'YLim', [0.5, length(sprois) + 0.5]);
set(gca, 'XTick', [], 'YTick', []);
xs = get(gca, 'XLim');
ys = get(gca, 'YLim');
colorbar;
axis square;

text(xs(1) + 0.6 * range(xs), ys(1) + 0.05 * range(ys), ...
     ['chi^2 test: p = ', sprintf('%.3f', p_chi2)]);

% --- Labels for columns --- %
ht_cols = nan(1, numel(sprois));
ht_cols_top = nan(1, numel(sprois));
for k1 = 1 : numel(sprois)
    ht_cols(k1) = text(k1, numel(sprois) + verticalPadding, ...
                       strrep(sprois{k1}, [hemi, '_'], ''), ...
                       'FontSize', 12);
    ht_cols_top(k1) = text(k1, 0, ...
                       strrep(sprois{k1}, [hemi, '_'], ''), ...
                       'FontSize', 12);
    set(ht_cols(k1), 'rotation', 90);
    set(ht_cols_top(k1), 'rotation', 90);
end

% --- Labels for rows --- %
ht_rows = nan(1, numel(sprois));
for k1 = 1 : numel(sprois)
    ht_rows(k1) = text(-horizontalPadding, k1, ...
                       strrep(sprois{k1}, [hemi, '_'], ''), ...
                       'FontSize', 12);
end


return