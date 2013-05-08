function imOut = fill_aparc12_roi_fig(imIn, imDiv, imText, hemi, fillROIs, ...
                                      varargin)
%% 
                                   
%% Process additional input arguments
bDEBUG = ~isempty(fsic(varargin, 'DEBUG'));

if ~isempty(fsic(varargin, 'fillClrs'))
    fillClrs = varargin{fsic(varargin, 'fillClrs') + 1};
else
    fillClrs = repmat({[1, 1, 0]}, 1, numel(fillROIs));
end


if ~iscell(fillClrs)
    error('fillClrs must be a cell structure');
end

if length(fillClrs) ~= length(fillROIs)
    error('Lengths of fillROIs and fillClrs do not match');
end

clrNameROIs = {};
clrNameROIClrs = {};
if ~isempty(fsic(varargin, 'clrNameROIs'))
    clrNameROIs = varargin{fsic(varargin, 'clrNameROIs') + 1};
    clrNameROIClrs = varargin{fsic(varargin, 'clrNameROIs') + 2};
end

drawBndROIs = {};
drawBndClrs = {};
if ~isempty(fsic(varargin, 'drawBndROIs'))
    drawBndROIs = varargin{fsic(varargin, 'drawBndROIs') + 1};
    drawBndClrs = varargin{fsic(varargin, 'drawBndROIs') + 2};
end

hashROIs = {};
hashROISigns = {};
if ~isempty(fsic(varargin, 'hashROIs'))
    hashROIs = varargin{fsic(varargin, 'hashROIs') + 1};
    hashROISigns = varargin{fsic(varargin, 'hashROIs') + 2};
end

if length(clrNameROIs) ~= length(clrNameROIClrs)
    error('Length mismatch between clrNameROIs and clrNameROIClrs');
end

%%
if ~isempty(hashROIs)
    hashBG = zeros(size(imIn, 1), size(imIn, 2), 1);
    
    hashSpc = 10;
    hashW = 2;
    for i1 = 1 : hashSpc : size(hashBG, 1)
        idx1 = i1;
        idx2 = 1;
        while (idx1 <= size(hashBG, 1)) && (idx2 <= size(hashBG, 2))
            hashBG(idx1, idx2 : idx2 + hashW) = 1;
            
            idx1 = idx1 + 1;
            idx2 = idx2 + 1;
        end
    end
    
    for i1 = 1 : hashSpc : size(hashBG, 2)
        idx1 = 1;
        idx2 = i1;
        while (idx1 <= size(hashBG, 1)) && (idx2 <= size(hashBG, 2))
            hashBG(idx1, idx2 : idx2 + hashW) = 1;
            
            idx1 = idx1 + 1;
            idx2 = idx2 + 1;
        end
    end
    
    hashBG = hashBG(1 : size(imIn, 1), 1 : size(imIn, 2));
    
    hashBGs{1} = hashBG;
    hashBGs{2} = fliplr(hashBG);
end

%%
imOut = imIn;

if length(size(imIn)) == 3
    imIn = imIn(:, :, 2);
end
if length(size(imDiv)) == 3
    imDiv = imDiv(:, :, 2);
end
if length(size(imText)) == 3
    imText = imText(:, :, 2);
end
imIn = 1 - double(imIn) / 255;
imDiv = 1 - double(imDiv) / 255;
imText = 1 - double(imText) / 255;

imto = imText - imIn;

if ~isequal(size(imIn), size(imDiv)) || ~isequal(size(imIn), size(imText))
    error('Size mismatch between the three input images');
end

imRegions = zeros(size(imDiv, 1), size(imDiv, 2), numel(fillROIs));
imRegionsClr = cell(1, numel(fillROIs));
bnds = cell(1, numel(fillROIs));

coords = [];
for i1 = 1 : numel(fillROIs)
    roiName = sprintf('%s.%s', hemi, fillROIs{i1});
    roiCoord = get_aparc12_roi_fig_coords(roiName, 1);
    t_imRegions = zeros(size(imDiv, 1), size(imDiv, 2), length(roiCoord));
    
    for i2 = 1 : length(roiCoord)
        if ~iscell(roiCoord)
            pause(0);
        end
        if iscell(roiCoord) && ~isempty(roiCoord)
            t_imRegions(:, :, i2) = regiongrowing(imDiv, roiCoord{i2}(2), roiCoord{i2}(1));
            coords = [coords; [roiCoord{i2}(2), roiCoord{i2}(1)]];
        else
            fprintf(2, 'WARNING: Skipping ROI %s, due to missing ROI coordinates.\n', ...
                    roiName)
        end
    end
    
    
    % -- Hash ROIs (optional) -- %
    if ~isempty(fsic(hashROIs, fillROIs{i1}))
        idxh = fsic(hashROIs, fillROIs{i1});
        if hashROISigns(idxh) < 0
            t_hashBG = hashBGs{1};
        else
            t_hashBG = hashBGs{2};
        end
        
        
        for j1 = 1 : size(t_imRegions, 3)
            t_im = t_imRegions(:, :, j1);
%         assert(isequal(size(t_imRegions), size(hashBG)));
            idx_h = find(t_im & t_hashBG);
        
            t_im(idx_h) = 0;
            
            t_imRegions(:, :, j1) = t_im;
        end
    end
    
    imRegions(:, :, i1) = sum(t_imRegions, 3);
    
    if ~isempty(fsic(drawBndROIs, fillROIs{i1}))
        imDil = imdilate(imRegions(:, :, i1), strel('disk', 3));
        bnd = bwboundaries(imDil);
        
        bnds{i1} = bnd{1};
    end
    
    % -- Furnish color information -- %
    imRegionsClr{i1} = zeros(size(imDiv, 1), size(imDiv, 2), 3);
    
    for k1 = 1 : 3
        t_clr_img = zeros(size(imDiv, 1), size(imDiv, 2));
        t_clr_img(find(imRegions(:, :, i1))) = fillClrs{i1}(k1);
        imRegionsClr{i1}(:, :, k1) = t_clr_img;
    end
end
imRegions = double(sum(imRegions, 3));


imOutB = imIn;
imOutB(find(imRegions)) = imRegions(find(imRegions));
imOut = repmat(imIn, [1, 1, 3]);

for i1 = 1 : numel(fillROIs)
    idxFill1 = find(imRegionsClr{i1}(:, :, 1));
    idxFill2 = find(imRegionsClr{i1}(:, :, 2));
    idxFill3 = find(imRegionsClr{i1}(:, :, 3));
    idxFill = union(idxFill1, idxFill2);
    idxFill = union(idxFill, idxFill3);
    
    for k1 = 1 : 3
        t_img = imOut(:, :, k1);
        t_clr_img = imRegionsClr{i1}(:, :, k1);
        t_img(idxFill) = 1.0 - t_clr_img(idxFill);
        
        imOut(:, :, k1) = t_img;
    end
end

imOut = 1.0 - imOut;

% imOut(:, :, 3) = imOutB;
% imOut = 1 - imOut;

% DEBUG;
if bDEBUG
    for i1 = 1 : size(coords, 1)
        imOut(coords(i1, 1), coords(i1, 2), :) = 0;
    end
end



%% (Optional): color the ROI names 
% if ~isempty(clrNameROIs)
imRegText = cell(1, length(clrNameROIs));

for i1 = 1 : numel(clrNameROIs)
    roiName = sprintf('%s.%s', hemi, clrNameROIs{i1});
    roiCoord = get_aparc12_roi_fig_coords(roiName, 1);
    t_imRegions = zeros(size(imDiv, 1), size(imDiv, 2), length(roiCoord));

    for i2 = 1 : length(roiCoord)
        if ~iscell(roiCoord)
            pause(0);
        end
        if iscell(roiCoord) && ~isempty(roiCoord)
            t_imRegions(:, :, i2) = regiongrowing(imDiv, roiCoord{i2}(2), roiCoord{i2}(1));
            coords = [coords; [roiCoord{i2}(2), roiCoord{i2}(1)]];
        else
            fprintf(2, 'WARNING: Skipping ROI %s, due to missing ROI coordinates.\n', ...
                    roiName)
        end
    end
    imRegText{i1} = t_imRegions;

end
    
    
idxText = find(imto > 0.5);
imtoClr_r = imto;
imtoClr_g = imto;
imtoClr_b = imto;

for i1 = 1 : numel(clrNameROIs)
    idxFill = find(imRegText{i1});
%         
%         idxFill = union(idxFill1, idxFill2);
%         idxFill = union(idxFill, idxFill3);

    idxFill = intersect(idxFill, idxText);

    imtoClr_r(idxFill) = 1 - clrNameROIClrs{i1}(1);
    imtoClr_g(idxFill) = 1 - clrNameROIClrs{i1}(2);
    imtoClr_b(idxFill) = 1 - clrNameROIClrs{i1}(3);

end

imto = cat(3, imtoClr_r, imtoClr_g, imtoClr_b);
% end

%% Optional: draw highlighting boundaries 
for i1 = 1 : numel(drawBndROIs)
    idx = fsic(fillROIs, drawBndROIs{i1});
    
    imbnd = zeros(size(imto));
    for i2 = 1 : size(bnds{idx}, 1)
        imbnd(bnds{idx}(i2, 1), bnds{idx}(i2, 2), :) = 1.0 - drawBndClrs{i1};
    end
    imbnd = imdilate(imbnd, strel('disk', 1));
    
    imto = imto + imbnd;
end

%% 
% imOut = imOut + repmat((-(imText - imIn)), [1, 1, 3]);
% imOut1 = nan(size(imOut, 1), size(imOut, 2), 0);
imOut_r = imOut(:, :, 1);
imOut_g = imOut(:, :, 2);
imOut_b = imOut(:, :, 3);

imto_r = imto(:, :, 1);
imto_g = imto(:, :, 2);
imto_b = imto(:, :, 3);

idx = union(find(imto_r > 0.2), find(imto_g > 0.2));
idx = union(idx, find(imto_b > 0.2));

imOut_r(idx) = 1 - imto_r(idx);
imOut_g(idx) = 1 - imto_g(idx);
imOut_b(idx) = 1 - imto_b(idx);
    
imOut = cat(3, imOut_r, imOut_g, imOut_b);

% figure;
% imshow(imOut); 
return