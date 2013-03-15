function imOut = fill_aparc12_roi_fig(imIn, imDiv, imText, hemi, fillROIs, varargin)
bDEBUG = ~isempty(fsic(varargin, 'DEBUG'))

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

if ~isequal(size(imIn), size(imDiv)) || ~isequal(size(imIn), size(imText))
    error('Size mismatch between the three input images');
end

imRegions = nan(size(imDiv, 1), size(imDiv, 2), numel(fillROIs));
coords = [];
for i1 = 1 : numel(fillROIs)
    
    roiCoord = get_aparc12_roi_fig_coords(sprintf('%s.%s', hemi, fillROIs{i1}), 1);
    t_imRegions = nan(size(imDiv, 1), size(imDiv, 2), length(roiCoord));
    
    for i2 = 1 : length(roiCoord)
        t_imRegions(:, :, i2) = regiongrowing(imDiv, roiCoord{i2}(2), roiCoord{i2}(1));
        coords = [coords; [roiCoord{i2}(2), roiCoord{i2}(1)]];
    end
    imRegions(:, :, i1) = sum(t_imRegions, 3);
end
imRegions = double(sum(imRegions, 3));

imOutB = imIn;
imOutB(find(imRegions)) = imRegions(find(imRegions));
imOut = repmat(imIn, [1, 1, 3]);
imOut(:, :, 3) = imOutB;
imOut = 1 - imOut;

% DEBUG;
if bDEBUG
    for i1 = 1 : size(coords, 1)
        imOut(coords(i1, 1), coords(i1, 2), :) = 0;
    end
end

imOut = imOut + repmat((-(imText - imIn)), [1, 1, 3]);

% figure;
% imshow(imOut); 
return