%% Make the plots for the human tsne coronal view (PCA sweep)
clear variables
close all

% Load additional color data
ids = [2 3 5 10 20 100];
markersize = 30;
figsubnums = 'abcdefghijkl';
idnames = {'9861','10021','12876','14380','15496','15697'};

%% Figure a-f
% t-SNE mapping colored by region with 10 components
for idnr = 1:length(ids)
% Load data
load(['results_human_average_pca_sweep/MappedHumanGenesAverage2Did' sprintf('%d',ids(idnr))],'mappedX','color_RGB');

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
set(gcf,'PaperPositionMode','auto')
saveas(1,['Figure8' figsubnums(idnr) '.png']);
end

%% Figure g-l
% t-SNE mapping colored by region with 10 components
for idnr = 1:length(ids)
% Load data
load(['results_human_average_pca_sweep/MappedHumanGenesAverage2Did' sprintf('%d',ids(idnr))],'mappedX','brain_id');

cm = hsv(max(brain_id));

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
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',cm(brain_id(snr),:),'Markersize',markersize);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
set(gcf,'PaperPositionMode','auto')
saveas(1,['Figure8' figsubnums(idnr+6) '.png']);
end

%% Figure legend
fs = 20;
figure(13)
clf
hold on
for idnr = 1:length(idnames)
    h = plot(0,idnr,'.','Color',cm(idnr,:),'MarkerSize',40);
    set(h,'visible','off')
end
h = legend(idnames,'Orientation','horizontal','FontSize',fs);
axis off
set(h,'Units','Pixel');
p = get(h,'Position');
set(h,'Position',[1 1 p(3) p(4)]);
pf = get(gcf,'Position');
set(gcf,'Position',[pf(1) pf(2) p(3) p(4)]);
set(h,'XColor',[1 1 1],'YColor',[1 1 1])
set(gcf,'PaperpositionMode','auto')
saveas(gcf,'Figure8_braincolors.eps','epsc');

