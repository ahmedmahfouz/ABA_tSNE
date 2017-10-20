%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');

%% Load data (currently for a single brain)
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
expressionlevels = load([datadir 'expressionMatrix']);
ngenes = size(expressionlevels.expressionMatrix,2);

%% Load gene selection
selection = load([datadir 'highConfGenesInd']);

% Load gene coordinates
coordsdata = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],coordsdata.voxel.brain_idx);
coords = [cx cy cz];

%% Perform t-SNE
X = expressionlevels.expressionMatrix;

selvec = randperm(size(X,1));
X = X(selvec,selection.highConfGenes);
coords = coords(selvec,:);
perplexity = 30;
theta = 0.7;
initial_dims = min(300,size(X,2));
mappedX = fast_tsne(X,3,initial_dims,perplexity,theta);

%% Save mapped data
save('MappedMouseGenesSelectionRand','X','mappedX','coords',...
    'initial_dims','perplexity','theta','selvec');

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

saveas(1,'MappedMouseGeneCoordsSelectionRand');

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

saveas(2,'MappedMouseGenesSelectionRand');

%% Image volume
isz = [67 41 58];
imgr = zeros(isz);
imgg = zeros(isz);
imgb = zeros(isz);

ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
imgr(ind) = (mappedX(:,1)-minx(1))./xrange(1);
imgg(ind) = (mappedX(:,2)-minx(2))./xrange(2);
imgb(ind) = (mappedX(:,3)-minx(3))./xrange(3);

save('MappedMouseGeneVolumeSelectionRand','imgr','imgg','imgb');
