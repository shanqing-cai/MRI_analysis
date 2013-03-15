function coord = get_aparc12_roi_fig_coords(roi, bFull)
dat.rh = {'OC',   {[27, 239], [27, 209], [52, 209]};
          'AG',   {[77, 175]};
          'MTO',  {[84, 248]};
          'ITO',  {[90, 289]};
          'SPL',  {[144, 92]};
          'pSMg', {[120, 157]};
          'aSMg', {[180, 167]};
          'dSC',  {[189, 85], [746, 92]};
          'vSC',  {[213, 171]};
          'dMC',  {[227, 73], [704, 79]};
          'midMC',{[236, 110]};
          'vMC',  {[239, 172]};
          'pdPMC',{[245, 71], [685, 65]};
          'midPMC', {[265, 113]};
          'vPMC', {[265, 171]};
          'mdPMC',{[268, 70], [275, 53]};
          'adPMC',{[301, 59]};
          'pMFg', {[298, 109], [320, 135]};
          'pIFs', {[300, 144], [290, 144]};
          'dIFo', {[302, 182]};
          'vIFo', {[292, 224]};
          'SFg',  {[344, 85], [379, 110], [526, 113]};
          'aMFg', {[350, 125], [338, 110]};
          'aIFs', {[341, 151]};
          'dIFt', {[352, 196]};
          'vIFt', {[346, 245]};
          'FOC',  {[346, 258], [345, 273]};
          'FP',   {[416, 199], [495, 185]}; %  [502, 105]
          'pSTg', {[192, 221]};
          'pdSTs',{[192, 247]};
          'pvSTs',{[192, 265]};
          'pMTg', {[192, 289]};
          'pITg', {[192, 310]};
          'aSTg', {[280, 248]};
          'adSTs',{[280, 272]};
          'avSTs',{[280, 293]};
          'aMTg', {[280, 311]};
          'aITg', {[280, 328]};
          'TP',   {[317, 308], [317, 258], [573, 310]};
          'PO',   {[419, 342]};
          'PT',   {[398, 420]};
          'pCO',  {[453, 353]};
          'pINS', {[457, 398]};
          'H',    {[440, 421]};
          'aCO',  {[475, 349]};
          'aINS', {[493, 382]};
          'PP',   {[502, 442]};
          'pFO',  {[504, 355]};
          'aFO',  {[521, 368]};
          'SFg',  {[537, 114]};
          'preSMA', {[596, 99]};
          'SMA',  {[630, 73]};
          'PCN',  {[784, 124]};
          'aCG', {[553, 146]};
          'pCG', {[721, 136]};
          'FMC',  {[521, 236]};
          'SCC',  {[583, 252]};
          'pPH',  {[692, 285]};
          'aPH',  {[662, 285]};
          'PHg',  {[662, 285]};
          'LG',   {[782, 245]};
          'aTFg', {[638, 326]};
          'pTFg', {[710, 306]};
          'TOF',  {[795, 276]}};
      
if nargin == 1
    bFull = 0;
end
      
dat.lh = cell(size(dat.rh));
for i1 = 1 : size(dat.lh, 1)
    dat.lh{i1, 1} = dat.rh{i1, 1};
    dat.lh{i1, 2} = dat.rh{i1, 2};
    for i2 = 1 : length(dat.lh{i1, 2})
        dat.lh{i1, 2}{i2}(1) = 895 - dat.lh{i1, 2}{i2}(1);
        dat.lh{i1, 2}{i2}(2) = dat.lh{i1, 2}{i2}(2) - 10;
    end
end

hemi = roi(1 : 2);
roi = roi(4 : end);
flds = fields(dat);
if isempty(fsic(flds, hemi))
    error('Unrecognized hemisphere: %s', hemi);
end

idx = fsic(dat.(hemi)(:, 1), roi);
if isempty(idx)
    coord = [NaN, NaN];
else
    coord = dat.(hemi){idx, 2};
end

if ~bFull
    coord = coord{1};
end
return