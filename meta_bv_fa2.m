function meta_bv_fa2(tfbv, measName, varargin)
colors.PWS = 'r';
colors.PFS = 'b';

tfbv_PWS = tfbv.PWS;
tfbv_PFS = tfbv.PFS;

N = size(tfbv_PWS, 2);

figure('Color', 'w', 'Name', measName, 'NumberTitle', 'off');
grps = fields(colors);
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    plot(linspace(0, 1, N), mean(tfbv.(grp)), '-', 'Color', colors.(grp));
    hold on;
    plot(linspace(0, 1, N), ...
         mean(tfbv.(grp)) + std(tfbv.(grp)) / sqrt(size(tfbv.(grp), 1)), ...
         '--', 'Color', colors.(grp));
    plot(linspace(0, 1, N), ...
         mean(tfbv.(grp)) - std(tfbv.(grp)) / sqrt(size(tfbv.(grp), 1)), ...
         '--', 'Color', colors.(grp));
end
xlabel('Normalize position');
ylabel('By-voxel FA');
set(gca, 'YLim', [0, 1]);

title(measName);

% -- Optional: by-voxel correlation with some measure --
if ~isempty(fsic(varargin, 'corr'))
    lc_k1 = nan(1, N);
    lc_p = nan(1, N);
    lc_r2 = nan(1, N);
    covar = varargin{fsic(varargin, 'corr') + 1};
    if size(covar.(grps{1}), 1) < size(covar.(grps{2}), 2)
        covar.(grps{1}) = transpose(covar.(grps{1}));
        covar.(grps{2}) = transpose(covar.(grps{2}));
    end
    covar = [covar.(grps{1}); covar.(grps{2})];
    for i1 = 1 : N
        x = [tfbv.(grps{1})(:, i1); tfbv.(grps{2})(:, i1)];
        [k, r2, p] = lincorr(x, covar);
        lc_k1(i1) = k(2);
        lc_p(i1) = p;
        lc_r2(i1) = r2;
    end
end
return