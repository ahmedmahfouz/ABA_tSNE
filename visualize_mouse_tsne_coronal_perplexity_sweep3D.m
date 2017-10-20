%% Make the plots for the mouse tsne coronal view
clear variables
close all

% Load additional color data
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
colors = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],colors.voxel.brain_idx);
coords = [cx cy cz];
perp = '256';

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load(['results_mouse_genes_perplexity_sweep_th0p2/MappedMouseGenesPerplexitySweep3Dperp' perp '.mat'],'mappedX','mappedRGB');

% Plot figure
cm = jet(58/2+1);
nsamples = size(mappedX,1);
figure(1)
clf
hold on
for snr = 1:nsamples
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',cm(abs(coords(snr,3)-58/2)+1,:));
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(1,['MappedMousePerplexitySweepWidthMirrored3Dperp' perp '.png']);
saveas(1,['MappedMousePerplexitySweepWidthMirrored3Dperp' perp '.fig']);

%% Same plot with regions
load(['results_mouse_genes_perplexity_sweep_th0p2/MappedMouseGenesPerplexitySweep3Dperp' perp '.mat'],'mappedX','mappedRGB');

% Plot figure
nsamples = size(mappedX,1);
figure(2)
clf
hold on
colorRGB = zeros(nsamples,3);
for snr = 1:nsamples
    if ~isempty(colors.voxel.color_HEX{snr})
        colorRGB(snr,:) = rgbconv(colors.voxel.color_HEX{snr});
    else
        colorRGB(snr,:) = [0 0 0];
    end
    plot3(mappedX(snr,1),mappedX(snr,2),mappedX(snr,3),'.','Color',colorRGB(snr,:));
end

% Format the plot
axis off equal tight
view(3)
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(2,['MappedMousePerplexitySweepRegion3D2perp' perp '.png']);
saveas(2,['MappedMousePerplexitySweepRegion3Dperp' perp '.fig']);

%% Same plot with regions
load(['results_mouse_genes_perplexity_sweep_th0p2/MappedMouseGenesPerplexitySweep3Dperp' perp '.mat'],'mappedX','mappedRGB');

% Plot figure
nsamples = size(mappedX,1);
figure(3)
clf
hold on
cm = (mappedX-min(mappedX(:)))/(max(mappedX(:))-min(mappedX(:)));
for snr = 1:nsamples
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',cm(snr,:));
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(3,['MappedMousePerplexitySweepMappedColors3Dperp' perp '.png']);
saveas(3,['MappedMousePerplexitySweepMappedColors3Dperp' perp '.fig']);

%% Interpolated plot
caxis = linspace(min(mappedX(:)),max(mappedX(:)),200);

[mx,my,mz] = meshgrid(caxis,caxis,caxis);

dt = DelaunayTri(mappedX);
[id,dist] = nearestNeighbor(dt,[mx(:) my(:) mz(:)]);

%% Put colors in space
imgr = reshape(colorRGB(id,1),[1 1 1]*length(caxis));
imgg = reshape(colorRGB(id,2),[1 1 1]*length(caxis));
imgb = reshape(colorRGB(id,3),[1 1 1]*length(caxis));
dth = 1;
imgr(dist>dth) = 0;
imgg(dist>dth) = 0;
imgb(dist>dth) = 0;

%% Select and show slice
slice = 80;
cimg = cat(3,imgr(:,:,slice),imgg(:,:,slice),imgb(:,:,slice));
figure(3)
imagesc(cimg)
axis equal off
set(gcf,'Color','white')

%% Volume with 3D images
ind = colors.voxel.brain_idx;
imgr = zeros([67 41 58]);
imgg = zeros([67 41 58]);
imgb = zeros([67 41 58]);

cm = (mappedX-min(mappedX(:)))/(max(mappedX(:))-min(mappedX(:)));
imgr(ind) = cm(:,1);
imgg(ind) = cm(:,2);
imgb(ind) = cm(:,3);

%% Select and show slice
slice = 20;
cimg = cat(3,imgr(:,:,slice),imgg(:,:,slice),imgb(:,:,slice));
figure(4)
imagesc(cimg)
axis equal off
set(gcf,'Color','white')

save(['results_mouse_genes_perplexity_sweep_th0p2/MappedMousePerplexitySweepVolume3Dperp' perp],'imgr','imgg','imgb');

