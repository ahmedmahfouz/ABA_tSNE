%% Make the plots for the mouse tsne coronal view
clear variables
close all

% Load additional color data
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
colors = load([datadir 'voxels']);

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did10.mat','mappedX','mappedRGB');

% Plot figure
nsamples = size(mappedX,1);
figure(1)
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
saveas(1,'Figure1a.png');

%% Figure b
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did10.mat');

% Scale LAB
imgr = imgr*2/3;
imgg = imgg*2/3;
imgb = imgb*2/3;

all0 = imgr==0 & imgg==0 & imgb==0;
imgr(all0) = 1;
imgg(all0) = 1;
imgb(all0) = 1;

tslice = 30;
sslice = 20;
cslice = 40;
figure(2)
timg = cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:)));
simg = cat(3,squeeze(imgr(:,:,sslice)),squeeze(imgg(:,:,sslice)),squeeze(imgb(:,:,sslice)));
cimg = cat(3,squeeze(imgr(cslice,:,:)),squeeze(imgg(cslice,:,:)),squeeze(imgb(cslice,:,:)));
img = cat(1,cat(2,cimg,flipdim(permute(simg,[2 1 3]),2)),cat(2,timg,ones(size(timg,1),size(simg,1),3)));

image(img)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(img,1)*10 size(img,2)*10],'Color','White');
saveas(2,'Figure1b.png');


%% Figure b sub
% t-SNE mapping colored by L*a*b
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did10.mat','mappedX','mappedRGB');

% Plot figure
nsamples = size(mappedX,1);
figure(3)
clf
hold on
for snr = 1:nsamples
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB(snr,:)*2/3);
end
% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(3,'Figure1bsub.png');

%% Figure d
% PCA mapping colored by region
load('results_mouse_highConf_coronal_SVD/SVDMouseGenesSelectionMeanSub.mat');

% Plot figure
nsamples = size(mappedX,1);
figure(4)
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
saveas(4,'Figure1d.png');

%% Figure e
load('results_mouse_highConf_coronal_SVD/SVDMouseGeneVolumeSelectionMeanSub');
isz = size(imgr);
mappedab = [mappedX(:,1)/max(abs(mappedX(:,1)))*127 mappedX(:,2)/max(abs(mappedX(:,2)))*127];
tC = makecform('lab2srgb');
mappedRGB = applycform([50*ones(size(mappedab,1),1) mappedab],tC);
imgr2 = zeros(isz);
imgg2 = zeros(isz);
imgb2 = zeros(isz);
imgr2(colors.voxel.brain_idx) = mappedRGB(:,1);
imgg2(colors.voxel.brain_idx) = mappedRGB(:,2);
imgb2(colors.voxel.brain_idx) = mappedRGB(:,3);

all0 = imgr==0 & imgg==0 & imgb==0;
imgr2(all0) = 1;
imgg2(all0) = 1;
imgb2(all0) = 1;

tslice = 30;
sslice = 20;
cslice = 40;
figure(5)
timg = cat(3,squeeze(imgr2(:,tslice,:)),squeeze(imgg2(:,tslice,:)),squeeze(imgb2(:,tslice,:)));
simg = cat(3,squeeze(imgr2(:,:,sslice)),squeeze(imgg2(:,:,sslice)),squeeze(imgb2(:,:,sslice)));
cimg = cat(3,squeeze(imgr2(cslice,:,:)),squeeze(imgg2(cslice,:,:)),squeeze(imgb2(cslice,:,:)));
img = cat(1,cat(2,cimg,flipdim(permute(simg,[2 1 3]),2)),cat(2,timg,ones(size(timg,1),size(simg,1),3)));

image(img)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(img,1)*10 size(img,2)*10],'Color','White');
saveas(5,'Figure1e.png');

%% Figure e sub
% t-SNE mapping colored by L*a*b
% Load data
load('results_mouse_highConf_coronal_SVD/SVDMouseGenesSelectionMeanSub','mappedX');
mappedRGB = [imgr2(colors.voxel.brain_idx),imgg2(colors.voxel.brain_idx),imgb2(colors.voxel.brain_idx)];

% Plot figure
nsamples = size(mappedX,1);
figure(6)
clf
hold on
for snr = 1:nsamples
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB(snr,:)*2/3);
end
% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
saveas(6,'Figure1esub.png');
