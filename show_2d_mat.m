function [nSigs, sigConnections, sigVals] = ...
    show_2d_mat(sig_rs, sprois, hemi, figName, visParams, varargin)
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
sigVals = [];

%% --- Grid --- %%
xs = get(gca, 'XLim');
ys = get(gca, 'YLim');

gridClr = [0.8, 0.8, 0.8];
% -- Vertical -- %
for x0 = xs(1) : 1.0 : xs(2)
    plot([x0, x0], [x0, ys(2)], '-', 'Color', gridClr);
%     plot([x0, x0], ys, '-', 'Color', gridClr);
end

% -- Horizontal -- %
for y0 = ys(1) : 1.0 : ys(2) - 1.0
    plot([xs(1), y0], [y0, y0], '-', 'Color', gridClr);
%     plot(xs, [y0, y0], '-', 'Color', gridClr);
end

plot(xs, xs, '-', 'Color', [0.5, 0.5, 0.5]);


%%
nSigs = [0, 0]; % [Pos, Neg]
for k1 = 1 : numel(sprois)
    for k2 = 1 : numel(sprois)
        if abs(sig_rs(k2, k1)) > abs(log10(0.005))
%             text(k1 - cellShift, k2 - cellShift, '*', 'Color', 'w');
            draw_asterisk(k1, k2, 'w');
        elseif abs(sig_rs(k2, k1)) > abs(log10(0.01))
%             text(k1 - cellShift, k2 - cellShift, 'X', 'Color', 'w');
            draw_x(k1, k2, 'w');
        elseif abs(sig_rs(k2, k1)) > abs(log10(0.05))
%             text(k1 - cellShift, k2 - cellShift, 'O', 'Color', 'w'); 
            draw_diamond(k1, k2, 'w');
            
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
            sigVals(end + 1) = sig_rs(k2, k1);
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

if isempty(fsic(varargin, 'noShowChi2Test'))
    text(xs(1) + 0.6 * range(xs), ys(1) + 0.05 * range(ys), ...
          ['chi^2 test: p = ', sprintf('%.3f', p_chi2)]);
end

% --- Labels for columns --- %
ht_cols = nan(1, numel(sprois));
ht_cols_top = nan(1, numel(sprois));
horizontalAdjust = -0.10;
verticalAdjust = 0.75;
for k1 = 1 : numel(sprois)
    ht_cols(k1) = text(k1, ...
                       numel(sprois) + verticalAdjust, ...
                       strrep(sprois{k1}, [hemi, '_'], ''), ...
                       'FontSize', 12);
    ht_cols_top(k1) = text(k1 + horizontalAdjust, 0, ...
                       strrep(sprois{k1}, [hemi, '_'], ''), ...
                       'FontSize', 12);
    set(ht_cols(k1), 'rotation', -90);
    set(ht_cols_top(k1), 'rotation', 90);
end

% --- Labels for rows --- %
ht_rows = nan(1, numel(sprois));
for k1 = 1 : numel(sprois)
    ht_rows(k1) = text(-horizontalPadding, k1, ...
                       strrep(sprois{k1}, [hemi, '_'], ''), ...
                       'FontSize', 12);
end

%% --- Draw legend --- %%
lgdX = 19;
lgdY0 = 3;
draw_diamond(lgdX, lgdY0, 'k');
text(lgdX + 1.5, lgdY0, 'p < 0.05');

draw_x(lgdX, lgdY0 + 2, 'k');
text(lgdX + 1.5, lgdY0 + 2, 'p < 0.01');

draw_asterisk(lgdX, lgdY0 + 4, 'k');
text(lgdX + 1.5, lgdY0 + 4, 'p < 0.005');

rectangle('Position', [lgdX - 1, lgdY0 - 1, 8, 6], ...          
          'EdgeColor', 'k', 'FaceColor', 'none');
      

%% --- Color bar labels --- %%
text(30, -1.5, 'Sig-value', 'Color', 'k');

if ~isempty(fsic(varargin, 'colorBarBGC'))
    text(30, 0, 'PWS>PFS', 'Color', 'r');
    text(30, 29, 'PWS<PFS', 'Color', 'b');
elseif ~isempty(fsic(varargin, 'colorBarCorr'))
    text(30, 0, 'Positive', 'Color', 'r');
    text(30, 29, 'Negative', 'Color', 'b');
end

return

%% Sub functions
function draw_diamond(k1, k2, clr)
    plot(k1 + [0, 0.5], k2 + [0.5, 0], '-', 'Color', clr);
    plot(k1 + [-0.5, 0], k2 + [0, 0.5], '-', 'Color', clr);
    plot(k1 + [0, 0.5], k2 + [-0.5, 0], '-', 'Color', clr);
    plot(k1 + [-0.5, 0], k2 + [0, -0.5], '-', 'Color', clr);
return

function draw_x(k1, k2, clr)
    plot(k1 + [-0.5, 0.5], k2 + [-0.5, 0.5], '-', 'Color', clr);
    plot(k1 + [0.5, -0.5], k2 + [-0.5, 0.5], '-', 'Color', clr);
return

function draw_asterisk(k1, k2, clr)
    plot(k1 + [-0.5, 0.5], k2 + [-0.5, 0.5], '-', 'Color', clr);
    plot(k1 + [0.5, -0.5], k2 + [-0.5, 0.5], '-', 'Color', clr);
    plot(k1 + [0, 0], k2 + [-0.5, 0.5], '-', 'Color', clr);
    plot(k1 + [-0.5, 0.5], k2 + [0, 0], '-', 'Color', clr);
return
