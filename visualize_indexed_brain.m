% Visualization script for indexed brain
clear
close all

%% Load data
id = '9861';
datadir = [LUMCDATADIR 'tSNE_ABA\normalized\'];
load([datadir 'gene_indices_full_brain_' id]);
load([datadir 'results_human_average_cell_selection\' 'MappedHumanGenesAverage2Dneurons' id 'id247']);

%% Give color to each label (each sample)
nsamples = size(mappedX,1);
maxval = max(abs(mappedX(:)));
Lconst = 50;

labimga = zeros(size(indimg));
labimgb = zeros(size(indimg));
labimgL = zeros(size(indimg));

for snr = 1:nsamples
    mappedab = mappedX(snr,:);
    labimga(indimg==snr) = mappedab(1)*127/maxval;
    labimgb(indimg==snr) = mappedab(2)*127/maxval;
    labimgL(indimg==snr) = Lconst;
end

%% Transform ab colorspace to Lab
tC = makecform('lab2srgb');
rgb = applycform([labimgL(:) labimga(:) labimgb(:)],tC);
rgbimgr = reshape(rgb(:,1),size(indimg));
rgbimgg = reshape(rgb(:,2),size(indimg));
rgbimgb = reshape(rgb(:,3),size(indimg));

%% Visualize slice
% Transversal slice
tsel = 100;
tslice = zeros([size(indimg,1) size(indimg,2) 3]);
tslice(:,:,1) = squeeze(rgbimgr(:,:,tsel));
tslice(:,:,2) = squeeze(rgbimgg(:,:,tsel));
tslice(:,:,3) = squeeze(rgbimgb(:,:,tsel));
figure(1)
clf
pos0 = get(gcf,'Position');
imshow(permute(tslice,[2 1 3]));
pos = [pos0(1) pos0(2) size(tslice,1) size(tslice,2)];
set(gcf,'Position',pos);
set(gca,'Units','normalized','Position',[0 0 1 1]);

% Coronal slice
csel = 80;
cslice = zeros([size(indimg,1) size(indimg,3) 3]);
cslice(:,:,1) = squeeze(rgbimgr(:,csel,:));
cslice(:,:,2) = squeeze(rgbimgg(:,csel,:));
cslice(:,:,3) = squeeze(rgbimgb(:,csel,:));
figure(2)
clf
pos0 = get(gcf,'Position');
imshow(flipdim(permute(cslice,[2 1 3]),1));
pos = [pos0(1) pos0(2) size(cslice,1) size(cslice,2)];
set(gcf,'Position',pos);
set(gca,'Units','normalized','Position',[0 0 1 1]);

% Sagittal slice
ssel = 80;
sslice = zeros([size(indimg,2) size(indimg,3) 3]);
sslice(:,:,1) = squeeze(rgbimgr(ssel,:,:));
sslice(:,:,2) = squeeze(rgbimgg(ssel,:,:));
sslice(:,:,3) = squeeze(rgbimgb(ssel,:,:));
figure(3)
clf
pos0 = get(gcf,'Position');
imshow(flipdim(permute(sslice,[2 1 3]),1));
pos = [pos0(1) pos0(2) size(sslice,1) size(sslice,2)];
set(gcf,'Position',pos);
set(gca,'Units','normalized','Position',[0 0 1 1]);
