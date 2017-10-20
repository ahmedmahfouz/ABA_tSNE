%% Make the plots for the mouse tsne coronal view (PCA sweep)
clear variables
close all

% Load additional color data
ids = {'9861','10021','12876','14380','15496','15697'};
markersize = 30;
figsubnums = 'abcdef';

%% Figure a-f
% t-SNE mapping colored by region with 10 components
for idnr = 1:length(ids)
% Load data
load(['results_human_average_pca_sweep/MappedHumanGenesAverage2D' ids{idnr} 'id300'],'mappedX','color_RGB');

% Plot figure
nsamples = size(mappedX,1);
figure(1)
clf
hold on
for snr = 1:nsamples
%     if ~isempty(colors.voxel.color_HEX{snr})
%         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
%     else
%         colorRGB = [0 0 0];
%     end
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',color_RGB(snr,:),'Markersize',markersize);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
% saveas(1,['Figure5' figsubnums(idnr) '.png']);
end
