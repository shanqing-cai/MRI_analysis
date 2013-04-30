function [t_bgc, p_bgc, r_lc, p_lc] = ...
    meta_bgComp_linCorr(meas, measName, behav, behavName, p_thresh, ...
                        sprois, varargin)
% === Compare strength and correlation with behav, node by node === %
nrois = length(sprois);

t_bgc = nan(nrois, 1);
p_bgc = nan(nrois, 1);

r_lc = nan(nrois, 1);
p_lc = nan(nrois, 1);

testName = 'ttest2';
if ~isempty(fsic(varargin, 'testName'))
    testName = varargin{fsic(varargin, 'testName') + 1};
end

fprintf(1, '=== Significant differences in %s: p_thresh (unc.) = %f ===\n', ...
        behavName, p_thresh);
for i1 = 1 : nrois
    if isequal(testName, 'ttest2')
        [tt_h, tt_p, ~, tt_stats] = ...
            ttest2(meas.PWS(i1, :), meas.PFS(i1, :));
        
        t_bgc(i1) = tt_stats.tstat;
        p_bgc(i1) = tt_p;
    elseif isequal(testName, 'ranksum')
        [rs_p, ~, rs_stats] = ranksum(meas.PWS(i1, :), meas.PFS(i1, :));
        
        t_bgc(i1) = rs_stats.zval;
        p_bgc(i1) = rs_p;
    end
       
    [t_lc_k, t_lc_r2, t_lc_p] = lincorr(meas.PWS(i1, :), behav);
    r_lc(i1) = sqrt(t_lc_r2) * sign(t_lc_k(2));
    p_lc(i1) = t_lc_p;
    
    if p_bgc(i1) < p_thresh
        if t_bgc(i1) < 0
            diffDirStr = 'PWS < PFS';
        else
            diffDirStr = 'PWS > PFS';
        end
        
        if isequal(testName, 'ttest2')
            fprintf(1, '%s: t = %.3f (%s); p = %.6f\n', ...
                    sprois{i1}, t_bgc(i1), diffDirStr, p_bgc(i1));
        elseif isequal(testName, 'ranksum')
            fprintf(1, '%s: z = %.3f (%s); p = %.6f\n', ...
                    sprois{i1}, t_bgc(i1), diffDirStr, p_bgc(i1));
        end
    end        
end
fprintf(1, '\n');

fprintf(1, '=== Significant correlations b/w %s and %s: p_thresh = %f ===\n', ...
        measName, behavName, p_thresh);
for i1 = 1 : nrois
    if p_lc(i1) < p_thresh
        fprintf(1, '%s: r = %.3f; p = %.6f\n', ...
                sprois{i1}, r_lc(i1), p_lc(i1));
    end
end
fprintf(1, '\n');

return