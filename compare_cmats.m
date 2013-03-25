function [p, sgn, varargout] = ...
    compare_cmats(a_cmat, bFold, testName, varargin)
%               p_SSI4_spr, rho_SSI4_spr, ...
%           p_EHcomp_spr, rho_EHcomp_spr, ...
%           p_rnSV_spr, rho_rnSV_spr] = ...
nrois = size(a_cmat.PWS, 1);

p = nan(nrois, nrois);
sgn = nan(nrois, nrois);

if nargin ==  6 || nargin == 7
    SSI4 = varargin{1};
    EH_comp = varargin{2};
    rnSV = varargin{3};
    
    rho_SSI4_spr = nan(nrois, nrois); 
    p_SSI4_spr = nan(nrois, nrois); % p-values from Spearman correlation with SSI4 (PWS only)

    rho_EHcomp_spr = nan(nrois, nrois);
    p_EHcomp_spr = nan(nrois, nrois);

    rho_rnSV_spr = nan(nrois, nrois);
    p_rnSV_spr = nan(nrois, nrois);
end


bRandPerm = ~isempty(fsic(varargin, '--randPerm'));
if bRandPerm
    idx_rp = randperm(size(a_cmat.PFS, 3) + size(a_cmat.PWS, 3));
end

bRandPermPWS = ~isempty(fsic(varargin, '--randPermPWS'));
if bRandPermPWS
    idx_rp_pws = randperm(size(a_cmat.PWS, 3));
end

for i1 = 1 : nrois
    for i2 = 1 : nrois
        if i1 == i2
            continue;
        end
        
        if bFold && (i2 > i1)
            continue;
        end
        
        v_PFS = squeeze(a_cmat.PFS(i1, i2, :));
        v_PWS = squeeze(a_cmat.PWS(i1, i2, :));
        
        if bRandPerm
            a_vs = [v_PFS; v_PWS];
            a_vs = a_vs(idx_rp);
            v_PFS = a_vs(1: length(v_PFS));
            v_PWS = a_vs(length(v_PFS) + 1 : end);
        end
        
        if bRandPermPWS
            v_PWS = v_PWS(idx_rp_pws);
        end
        
        if isequal(testName, 'ttest2')
            [~, p(i1, i2)] = ttest2(v_PFS, v_PWS);
        elseif isequal(testName, 'ranksum')
            [p(i1, i2), ~, stats] = ranksum(v_PFS, v_PWS);
        else
            error('Unrecognized test name: %s', testName);
        end
            
        if median(v_PWS) < median(v_PFS)
            sgn(i1, i2) = -1;
        else
            sgn(i1, i2) = 1;
        end
        
        if nargin == 6 || nargin == 7
            % --- Correlation with SSI4 --- %
            [rho_SSI4_spr(i1, i2), ~, p_SSI4_spr(i1, i2)] = ...
                spear(v_PWS, SSI4(:));

            % --- Correlation with EH_comp --- %
            v_2grp = [v_PWS; v_PFS];
            EHcomp_2grp = [EH_comp.PWS, EH_comp.PFS];
            [rho_EHcomp_spr(i1, i2), ~, p_EHcomp_spr(i1, i2)] = ...
                spear(v_2grp(:), EHcomp_2grp(:));

            % --- Correlation with rnSV --- %
            rnSV_2grp = [rnSV.PWS, rnSV.PFS];
            [rho_rnSV_spr(i1, i2), ~, p_rnSV_spr(i1, i2)] = ...
                spear(v_2grp(:), rnSV_2grp(:));
        end
    end
end

if nargin == 6 || nargin == 7
    varargout{1} = p_SSI4_spr;
    varargout{2} = rho_SSI4_spr;
    varargout{3} = p_EHcomp_spr;
    varargout{4} = rho_EHcomp_spr;
    varargout{5} = p_rnSV_spr;
    varargout{6} = rho_rnSV_spr;
end


return