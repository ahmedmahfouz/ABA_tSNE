%% Make the plots for the mouse tsne coronal view
clear variables
close all

% Load additional color data
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_sagittal/';
colors = load([datadir 'voxels']);

%% Figure a
% t-SNE mapping colored by region with 10 components
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Dsagittalid10.mat','mappedX','mappedRGB');

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
saveas(1,'Figure4a.png');

%% Figure b
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Dsagittalid10.mat');

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
saveas(2,'Figure4b.png');


%% Figure bsub
% t-SNE mapping colored by L*a*b
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Dsagittalid10.mat','mappedX','mappedRGB');

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
saveas(3,'Figure4bsub.png');

%% Figure c
% 
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGeneVolumeSelection2Dsagittalid100.mat');

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
figure(4)
timg = cat(3,squeeze(imgr(:,tslice,:)),squeeze(imgg(:,tslice,:)),squeeze(imgb(:,tslice,:)));
simg = cat(3,squeeze(imgr(:,:,sslice)),squeeze(imgg(:,:,sslice)),squeeze(imgb(:,:,sslice)));
cimg = cat(3,squeeze(imgr(cslice,:,:)),squeeze(imgg(cslice,:,:)),squeeze(imgb(cslice,:,:)));
img = cat(1,cat(2,cimg,flipdim(permute(simg,[2 1 3]),2)),cat(2,timg,ones(size(timg,1),size(simg,1),3)));

image(img)
axis off tight equal
set(gca,'Units','normalize','Position',[0 0 1 1]);
set(gcf,'Units','pixel','Position',[50 50 size(img,1)*10 size(img,2)*10],'Color','White');
saveas(4,'Figure4c.png');

%% Figure csub
% t-SNE mapping colored by L*a*b
% Load data
load('results_mouse_genes_selection_pca_sweep/MappedMouseGenesSelection2Dsagittalid100.mat','mappedX','mappedRGB');

% Plot figure
nsamples = size(mappedX,1);
figure(5)
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
saveas(5,'Figure4csub.png');
