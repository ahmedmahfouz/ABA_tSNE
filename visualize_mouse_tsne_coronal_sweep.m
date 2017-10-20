%% Make the plots for the mouse tsne coronal view (PCA sweep)
clear variables
close all

% Load additional color data
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
colors = load([datadir 'voxels']);

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did2.mat','mappedX','mappedRGB');

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
saveas(1,'Figure2a.png');

%% Figure b
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did3.mat','mappedX','mappedRGB');

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
saveas(1,'Figure2b.png');

%% Figure c
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did5.mat','mappedX','mappedRGB');

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
saveas(1,'Figure2c.png');

%% Figure d
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
saveas(1,'Figure2d.png');

%% Figure e
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did20.mat','mappedX','mappedRGB');

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
saveas(1,'Figure2e.png');

%% Figure f
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did100.mat','mappedX','mappedRGB');

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
saveas(1,'Figure2f.png');

%% Figure g
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did2.mat');

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
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2g.png');

%% Figure h
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did3.mat');

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
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2h.png');

%% Figure i
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did5.mat');

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
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2i.png');

%% Figure j
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
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2j.png');

%% Figure k
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did20.mat');

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
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2k.png');

%% Figure l
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did100.mat');

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
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2l.png');

%% Figure m
% t-SNE mapping colored by slice
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Did100.mat','mappedX','mappedRGB','coords');

rng(0);
cm = rand(67,3)/2;

% Plot figure
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
saveas(1,'Figure2_slice_based_coloring_tsne.png');

%% Figure n
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Did100.mat');

% Scale LAB
imgr = imgr*2/3;
imgg = imgg*2/3;
imgb = imgb*2/3;

all0 = imgr==0 & imgg==0 & imgb==0;

for slnr = 1:max(coords(:,1))
    tmp0 = all0(slnr,:,:);
    tmp = imgr(slnr,:,:);
    tmp(~tmp0) = cm(slnr,1);
    imgr(slnr,:,:) = tmp;
    tmp = imgg(slnr,:,:);
    tmp(~tmp0) = cm(slnr,2);
    imgg(slnr,:,:) = tmp;
    tmp = imgb(slnr,:,:);
    tmp(~tmp0) = cm(slnr,3);
    imgb(slnr,:,:) = tmp;
end

imgr(all0) = 1;
imgg(all0) = 1;
imgb(all0) = 1;

tslice = 30;
sslice = 20;
cslice = 40;
figure(2)
timg = permute(cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:))),[2 1 3]);

imshow(timg)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(timg,2)*10 size(timg,1)*10],'Color','White');
saveas(2,'Figure2_slice_based_coloring_anatomical.png');