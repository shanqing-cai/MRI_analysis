function plot_elisa_speech_network()
hemi = 'lh';

roiFigs.lh = sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc_noText.tif', hemi);
roiFigsDiv.lh = sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc_noText_div.tif', hemi);
roiFigsText.lh = sprintf('/users/cais/STUT/figures/rois_%s_flat_SLaparc.tif', hemi);

snROIs = get_aparc12_cortical_rois('speech');

check_file(roiFigs.lh);
check_file(roiFigsDiv.lh);
check_file(roiFigsText.lh);

im = imread(roiFigs.lh);
imDiv = imread(roiFigsDiv.lh);
imText = imread(roiFigsText.lh);

imFilled = fill_aparc12_roi_fig(im, imDiv, imText, hemi, snROIs);

figFN = sprintf('/users/cais/STUT/figures/SpeechNetwork_elisa.tif');

hf = figure('Color', 'w');
imshow(imFilled);
saveas(hf, figFN, 'tif');
check_file(figFN);
fprintf('Saved image to %s\n', figFN);
return