%% Make the plots for the mouse tsne coronal view
clear variables
close all

% Load additional color data
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_sagittal/';
colors = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],colors.voxel.brain_idx);
coords = [cx cy cz];
perp = '256';

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load(['results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Dsagittalperp' perp '.mat'],'mappedX','mappedRGB');

% Plot figure
cm = jet(58/2+1);
cm = jet(67);
nsamples = size(mappedX,1);
figure(1)
clf
hold on
for snr = 1:nsamples
%     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',cm(abs(coords(snr,3)-58/2)+1,:));
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',cm(coords(snr,1),:));
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(1,['MappedMousePerplexitySweepSagittalSlicenumber2Dperp' perp '.png']);

%% Same plot with regions
load(['results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Dsagittalperp' perp '.mat'],'mappedX','mappedRGB');

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
saveas(2,['MappedMousePerplexitySweepSagittalRegion2Dperp' perp '.png']);

%% 3D volume
imgr = zeros([67 41 58]);
imgg = zeros([67 41 58]);
imgb = zeros([67 41 58]);

minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);
mappedab = mappedX*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([50*ones(size(mappedab,1),1) mappedab],tC);

ind = colors.voxel.brain_idx;
imgr(ind) = mappedRGB(:,1);
imgg(ind) = mappedRGB(:,2);
imgb(ind) = mappedRGB(:,3);

%% Show slice
slice = 25;
figure(3)
imshow(cat(3,squeeze(imgr(:,slice,:)),...
    squeeze(imgg(:,slice,:)),...
    squeeze(imgb(:,slice,:))));

save(['MappedMousePerplexitySweepSagittalVolume2Dperp' perp],'imgr','imgg','imgb');




