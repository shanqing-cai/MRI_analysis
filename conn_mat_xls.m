function conn_mat_xls(matFN, grp, xlsFN)
dat = load(matFN);

assert(isfield(dat.a_cmat, grp));

nROIs = length(dat.sprois);
nEntries = (1 + nROIs) * nROIs / 2;

A = cell(nEntries + 1, 3);
A(1, :) = {'ROI1', 'ROI2', 'TMN'};


row = 2;
mn_cmat = nanmean(dat.a_cmat.(grp), 3);
for i1 = 1 : nROIs
    for i2 = i1 : nROIs
        A{row, 1} = deblank(dat.sprois{i1});
        A{row, 2} = deblank(dat.sprois{i2});
        A{row, 3} = mn_cmat(i1, i2);
        
        row = row + 1;
    end
end


xlswrite(xlsFN, A, 1);
check_file(xlsFN);
return