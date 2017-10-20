%% Make the plots for the human tsne coronal view (PCA sweep)
clear variables
close all

% Load additional color data
markersize = 20;
figsubnums = 'abcdefghijkl';

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load('results_human_average_pca_sweep/MappedHumanGenesAverage2Did10','mappedX','color_RGB');

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
saveas(1,'Figure7a.png');

%% Figure b
% t-SNE mapping colored by region with 10 components
% Load data
load('results_human_average_group_SVD/SVDHumanGenesGroupAverageMeanSub','mappedX','color_RGB');

% Plot figure
nsamples = size(mappedX,1);
figure(2)
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
saveas(2,'Figure7b.png');

%% Figure c
% t-SNE mapping colored by region with 10 components
% Load data
load('results_human_average_group_MDS/MDSHumanGenesGroupAverageMeanSub','mappedX2','color_RGB');

% Plot figure
nsamples = size(mappedX2,1);
figure(3)
clf
hold on
for snr = 1:nsamples
%     if ~isempty(colors.voxel.color_HEX{snr})
%         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
%     else
%         colorRGB = [0 0 0];
%     end
    plot(mappedX2(snr,1),mappedX2(snr,2),'.','Color',color_RGB(snr,:),'Markersize',markersize);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
set(gcf,'PaperPositionMode','auto')
saveas(3,'Figure7c.png');

%% Figure d
% t-SNE mapping colored by region with 10 components
% Load data
load('results_human_average_pca_sweep/MappedHumanGenesAverage2Did10','mappedX','color_RGB','brain_id');
nbrains = max(brain_id);
cm = hsv(nbrains);

% Plot figure
figure(1)
clf
hold on
for bnr = 1:nbrains
%     if ~isempty(colors.voxel.color_HEX{snr})
%         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
%     else
%         colorRGB = [0 0 0];
%     end
    plot(mappedX(brain_id==bnr,1),mappedX(brain_id==bnr,2),'.','Color',cm(bnr,:),'Markersize',markersize);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
set(gcf,'PaperPositionMode','auto')
saveas(1,'Figure7d.png');

%% Figure e
% t-SNE mapping colored by region with 10 components
% Load data
load('results_human_average_group_SVD/SVDHumanGenesGroupAverageMeanSub','mappedX','color_RGB');

% Plot figure
figure(2)
clf
hold on
for bnr = 1:nbrains
%     if ~isempty(colors.voxel.color_HEX{snr})
%         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
%     else
%         colorRGB = [0 0 0];
%     end
    plot(mappedX(brain_id==bnr,1),mappedX(brain_id==bnr,2),'.','Color',cm(bnr,:),'Markersize',markersize);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
set(gcf,'PaperPositionMode','auto')
saveas(2,'Figure7e.png');

%% Figure f
% t-SNE mapping colored by region with 10 components
% Load data
load('results_human_average_group_MDS/MDSHumanGenesGroupAverageMeanSub','mappedX2','color_RGB');

% Plot figure
nsamples = size(mappedX2,1);
figure(3)
clf
hold on
for bnr = 1:nbrains
%     if ~isempty(colors.voxel.color_HEX{snr})
%         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
%     else
%         colorRGB = [0 0 0];
%     end
    plot(mappedX2(brain_id==bnr,1),mappedX2(brain_id==bnr,2),'.','Color',cm(bnr,:),'Markersize',markersize);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
set(gcf,'PaperPositionMode','auto')
saveas(3,'Figure7f.png');

