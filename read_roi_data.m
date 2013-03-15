function [roiNames, roiMat] = read_roi_data(roiDataTxt)                                             
    txt = textread(roiDataTxt, '%s', 'delimiter', '\n');

    roiNames = {};
    roiMat = [];

    for i1 = 1 : numel(txt)
        t_line = txt{i1};
        t_items = splitstring(t_line);
        t_vec = [];
        if str2num(deblank(t_items{2})) > 0
            roiNames{end + 1} = strrep(deblank(t_items{1}), ':', '');
            for i2 = 2 : length(t_items)
                t_vec(end + 1) = str2num(deblank(t_items{i2}));
            end

            roiMat(:, end + 1) = t_vec';
        else                                                                                 
            fprintf('Skipping ROI with all-zero data: %s\n', deblank(t_items{1}));                   
        end
    end
return
