%% Make the plots for the mouse tsne coronal view
clear variables
close all

% Load additional color data
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
colors = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],colors.voxel.brain_idx);
coords = [cx cy cz];
id = '100';

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load(['results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did' id '.mat'],'mappedX','mappedRGB');

% Plot figure
cm = rand(67,3)*0.8;
nsamples = size(mappedX,1);
figure(1)
clf
hold on
for snr = 1:nsamples
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',cm(coords(snr,1),:));
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(1,['MappedMouseGenesSelectionSlicenumber2Did' id '.png']);

%% Same plot with regions
load(['results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did' id '.mat'],'mappedX','mappedRGB');

% Plot figure
nsamples = size(mappedX,1);
figure(2)
clf
hold on
for snr = 1:nsamples
    if ~isempty(colors.voxel.color_HEX{snr})
        colorRGB = rgbconv(colors.voxel.color_HEX{snr});
    else
        colorRGB = [0 0 0];
    end
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',colorRGB);
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(2,['MappedMouseGenesSelectionRegion2Did' id '.png']);

