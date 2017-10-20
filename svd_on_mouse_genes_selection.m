%% Perform svd on mouse gene data
clear variables
close all

addpath('../bh-tsne');

%% Load data (currently for a single brain)
outdir = 'results_mouse_highConf_coronal_SVD/';
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
expressionlevels = load([datadir 'expressionMatrix']);
ngenes = size(expressionlevels.expressionMatrix,2);

%% Load gene selection
selection = load([datadir 'highConfGenesInd_coronal']);

% Load gene coordinates
coordsdata = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],coordsdata.voxel.brain_idx);
coords = [cx cy cz];

%% Perform t-SNE
X = expressionlevels.expressionMatrix;

ss = 1;
X = X(1:ss:end,selection.highConfGenes_coronal);
X = X-repmat(mean(X),[size(X,1) 1]);
coords = coords(1:ss:end,:);
covX = X' * X;
[M, lambda] = eig(covX);
[~, ind] = sort(diag(lambda), 'descend');
% if initial_dims > size(M, 2)
%     initial_dims = size(M, 2);
% end
M = M(:,ind);
mappedXc = X * M;
mappedX = mappedXc(:,1:3);

%% Save mapped data
save([outdir 'SVDMouseGenesSelectionMeanSub'],'X','mappedXc','mappedX','coords',...
    'M','lambda');

%% Visualize mapped data
% Colorspace
crange = max(coords)-min(coords);
mincoords = min(coords);
nsamples = size(X,1);

%
figure(1)
clf
hold on
if size(mappedX,2) == 2
    for snr = 1:nsamples
        c = (coords(snr,:)-mincoords)./crange;
        plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
    end
else
    for snr = 1:nsamples
        c = (coords(snr,:)-mincoords)./crange;
        plot3(mappedX(snr,1),mappedX(snr,2),mappedX(snr,3),'.','Color',c);
    end
    view(3)
    axis equal
end

saveas(1,[outdir 'SVDMouseGeneCoordsSelectionMeanSub']);

%% Visualize mapped data based on t-SNE coordinates
xrange = max(mappedX)-min(mappedX);
minx = min(mappedX);

figure(2)
clf
hold on
for snr = 1:nsamples
    c = (mappedX(snr,:)-minx)./xrange;
    cc = zeros(1,3);
    cc(1:length(c)) = c;
    plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',cc);
end
view(3)
axis equal

saveas(2,[outdir 'SVDMouseGenesSelectionMeanSub']);

%% Image volume
isz = [67 41 58];
imgr = zeros(isz);
imgg = zeros(isz);
imgb = zeros(isz);

imgr(coordsdata.voxel.brain_idx) = (mappedX(:,1)-minx(1))./xrange(1);
imgg(coordsdata.voxel.brain_idx) = (mappedX(:,2)-minx(2))./xrange(2);
imgb(coordsdata.voxel.brain_idx) = (mappedX(:,3)-minx(3))./xrange(3);

save([outdir 'SVDMouseGeneVolumeSelectionMeanSub'],'imgr','imgg','imgb');

%% Image with all components
cimg = zeros([isz size(mappedXc,2)]);
for cnr = 1:size(mappedXc,2)
    tmpimg = zeros(isz);
    tmpimg(coordsdata.voxel.brain_idx) = mappedXc(:,cnr);
    cimg(:,:,:,cnr) = tmpimg;
end

save([outdir 'SVDMouseGeneVolumeSelectionMeanSubAllcomponents'],'-v7.3','cimg');
