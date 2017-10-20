%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');

perplexity = 256;

%% Load data (currently for a single brain)
outdir = 'results_mouse_genes_selection_pca_sweep/';
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
expressionlevels = load([datadir 'expressionMatrix']);
ngenes = size(expressionlevels.expressionMatrix,2);
colors = load([datadir 'voxels']);

sagdatadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_sagittal/';
sagcolor = load([sagdatadir 'voxels']);

%% Load gene selection
selection = load([datadir 'highConfGenesInd_coronal']);

%% Check overlapping sample coordinates
imvol = [67 41 58];
cvox = false(imvol);
svox = false(imvol);
cvox(colors.voxel.brain_idx) = true;
svox(sagcolor.voxel.brain_idx) = true;
bvox = cvox&svox;
indimg = zeros(imvol);
indimg(colors.voxel.brain_idx) = 1:size(colors.voxel.brain_idx,1);
sampleselind = indimg(bvox);
[cx,cy,cz] = ind2sub(imvol,colors.voxel.brain_idx);
coords = [cx cy cz];

%% Perform t-SNE
X = expressionlevels.expressionMatrix;

% selvec = randperm(size(X,1));
sampleselvec = sampleselind;
geneselvec = selection.highConfGenes_coronal;
X = X(sampleselvec,geneselvec);

% Normalize
X = zscore(X);

coords = coords(sampleselvec,:);
theta = 0.2;
initial_dims = [];
mappedX = fast_tsne(X,2,initial_dims,perplexity,theta);

%
minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);
mappedab = mappedX*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);

%% Save mapped data
save([outdir sprintf('MappedMouseGenesSelection2DcoronalSelSagittalperp%d',perplexity)],'X','mappedX','coords',...
    'initial_dims','perplexity','theta','sampleselvec','geneselvec','mappedRGB');


% %% Visualize mapped data
% % Colorspace
% crange = max(coords)-min(coords);
% mincoords = min(coords);
% nsamples = size(X,1);
%
% %
% figure(1)
% clf
% hold on
% if size(mappedX,2) == 2
%     for snr = 1:nsamples
%         c = (coords(snr,:)-mincoords)./crange;
%         plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
%     end
% else
%     for snr = 1:nsamples
%         c = (coords(snr,:)-mincoords)./crange;
%         plot3(mappedX(snr,1),mappedX(snr,2),mappedX(snr,3),'.','Color',c);
%     end
%     view(3)
%     axis equal
% end
%
% saveas(1,'MappedMouseGeneCoordsSelection2D300');
% saveas(1,'MappedMouseGeneCoordsSelection2D300.png');

% %% Visualize mapped data based on t-SNE coordinates
% minx = min(mappedX);
% maxx = max(mappedX);
% maxval = max([abs(minx) abs(maxx)]);
% mappedab = mappedX*127/maxval;
% tC = makecform('lab2srgb');
% mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);
%
% figure(2)
% clf
% hold on
% for snr = 1:nsamples
%     plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',mappedRGB(snr,:));
% end
% view(3)
% axis equal
%
% saveas(2,'MappedMouseGenesSelection2D300');
% saveas(2,'MappedMouseGenesSelection2D300.png');

%% Visualize mapped data based on t-SNE coordinates
nsamples = size(mappedX,1);
figure(3)
clf
hold on
for snr = 1:nsamples
    if ~isempty(colors.voxel.color_HEX{sampleselvec(snr)})
        colorRGB = rgbconv(colors.voxel.color_HEX{sampleselvec(snr)});
    else
        colorRGB = [0 0 0];
    end
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',colorRGB);
end

saveas(3,[outdir sprintf('MappedMouseGeneSelection2DcoronalSelSagittalperp%d',perplexity)]);
saveas(3,[outdir sprintf('MappedMouseGeneSelection2DcoronalSelSagittalperp%d.png',perplexity)]);

%% Image volume
isz = [67 41 58];
imgr = zeros(isz);
imgg = zeros(isz);
imgb = zeros(isz);

ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
imgr(ind) = mappedRGB(:,1);
imgg(ind) = mappedRGB(:,2);
imgb(ind) = mappedRGB(:,3);

save([outdir sprintf('MappedMouseGeneVolumeSelection2DcoronalSelSagittalperp%d',perplexity)],'imgr','imgg','imgb');
