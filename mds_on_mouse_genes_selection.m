%% Perform svd on mouse gene data
clear variables
close all

addpath('../bh-tsne');

%% Load data (currently for a single brain)
outdir = 'results_mouse_highConf_coronal_MDS/';
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
expressionlevels = load([datadir 'expressionMatrix']);
ngenes = size(expressionlevels.expressionMatrix,2);
colors = load([datadir 'voxels']);

%% Load gene selection
selection = load([datadir 'highConfGenesInd_coronal']);

% Load gene coordinates
coordsdata = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],coordsdata.voxel.brain_idx);
coords = [cx cy cz];

%% Pre-process data
X = expressionlevels.expressionMatrix;

ss = 1;
X = X(:,selection.highConfGenes_coronal);
% Per gene
X = zscore(X,0,1);

%% Compute distance matrix
disp('Distance matrix');
nsamples = size(X,1);
D = zeros(nsamples,nsamples);
tic
parfor snr1 = 1:nsamples
    dist = zeros(1,nsamples);
    for snr2 = 1:nsamples
        dist(snr2) = sqrt(1-xcov(X(snr1,:),X(snr2,:),0,'coeff')^2);
    end
    D(snr1,:) = dist;
end
toc

D = abs(D);
D(eye(size(D))>0) = 0;

%% MDS
disp('2D MDS');
mappedX = mdscale(D,2);
% disp('3D MDS');
% mappedX3 = mdscale(D,3);

minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);
mappedab = mappedX*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);

%% Save mapped data
save([outdir 'MDSMouseGenesSelection'],'X','mappedX','coords',...
    'M','D','mappedRGB');

%% Visualize mapped data
nsamples = size(mappedX,1);

%
figure(1)
clf
hold on

for snr = 1:nsamples
    if ~isempty(colors.voxel.color_HEX{snr})
        colorRGB = rgbconv(colors.voxel.color_HEX{snr});
    else
        colorRGB = [0 0 0];
    end
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB);
end

saveas(1,[outdir 'MDSMouseGeneRegionsSelection']);

%% Visualize mapped data based on t-SNE coordinates
xrange = max(mappedX)-min(mappedX);
minx = min(mappedX);

figure(2)
clf
hold on
for snr = 1:nsamples
    plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',mappedRGB(snr,:));
end
view(3)
axis equal

saveas(2,[outdir 'MDSMouseGenesSelection']);

%% Image volume
isz = [67 41 58];
imgr = zeros(isz);
imgg = zeros(isz);
imgb = zeros(isz);

ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
imgr(ind) = mappedRGB(:,1);
imgg(ind) = mappedRGB(:,2);
imgb(ind) = mappedRGB(:,3);

save([outdir 'MDSMouseGeneVolumeSelection'],'imgr','imgg','imgb');
