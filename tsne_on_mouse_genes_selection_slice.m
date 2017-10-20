%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');

%% Load data (currently for a single brain)
disp('Load data');
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_coronal/';
expressionlevels = load([datadir 'expressionMatrix']);
ngenes = size(expressionlevels.expressionMatrix,2);

%% Load gene selection
disp('Load selection')
selection = load([datadir 'highConfGenesInd_coronal']);

% Load gene coordinates
coordsdata = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],coordsdata.voxel.brain_idx);
coords = [cx cy cz];

%% Perform t-SNE
disp('Perform t-SNE');
ndims = 2;
X = expressionlevels.expressionMatrix(:,selection.highConfGenes_coronal);
mappedX = zeros(size(X,1),ndims);
ss = 1;
% Slice direction is x, with spacing 1
for xnr = min(cx):ss:max(cx)
    selvec = cx==xnr;
    Xsel = X(selvec,:);
    Xsel = zscore(Xsel);
    X(selvec,:) = Xsel;
    perplexity = 30;
    theta = 0.5;
    initial_dims = min(100,size(Xsel,2));
    mappedX(selvec,:) = fast_tsne(Xsel,ndims,initial_dims,perplexity,theta);
end

%% Save mapped data
disp('Save data');
save('MappedMouseGenesSelectionSlice2D','X','mappedX','coords',...
    'initial_dims','perplexity','theta');

%% Visualize mapped data
% Colorspace
crange = max(coords)-min(coords);
mincoords = min(coords);
nsamples = size(X,1);

%
for xnr = min(cx):ss:max(cx)
    figure(1)
    clf
    hold on
    selvec = cx==xnr;
    selind = 1:size(mappedX,1);
    selind = selind(selvec);
    if size(mappedX,2) == 2
        for snr = selind
            c = (coords(snr,:)-mincoords)./crange;
            plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
        end
    else
        for snr = selind
            c = (coords(snr,:)-mincoords)./crange;
            plot3(mappedX(snr,1),mappedX(snr,2),mappedX(snr,3),'.','Color',c);
        end
        view(3)
        axis equal
    end
    
%     saveas(1,'MappedMouseGeneCoordsSelectionSlice2D');
    saveas(1,sprintf('MappedMouseGeneCoordsSelectionSlice2D_%d',xnr));
end

%% Visualize mapped data based on t-SNE coordinates
minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);
mappedab = mappedX*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);

for xnr = min(cx):ss:max(cx)
    figure(2)
    clf
    hold on
    selvec = cx==xnr;
    selind = 1:size(mappedX,1);
    selind = selind(selvec);
    for snr = selind
        plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',mappedRGB(snr,:));
    end
    view(3)
    axis equal
    
    % saveas(2,'MappedMouseGenesSelectionSlice2D');
    saveas(2,sprintf('MappedMouseGenesSelectionSlice2D_%d',xnr));
end

%% Visualize mapped data based on t-SNE coordinates
minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);
mappedab = mappedX*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);

for xnr = min(cx):ss:max(cx)
    figure(3)
    clf
    hold on
    selvec = cx==xnr;
    selind = 1:size(mappedX,1);
    selind = selind(selvec);
    for snr = selind
        plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB(snr,:));
    end
    
    % saveas(3,'MappedMouseGenesSelectionLabMapSlice2D');
    saveas(3,sprintf('MappedMouseGenesSelectionLabMapSlice2D_%d.png',xnr));
end

%% Image volume
isz = [67 41 58];
imgr = zeros(isz);
imgg = zeros(isz);
imgb = zeros(isz);

ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
imgr(ind) = mappedRGB(:,1);
imgg(ind) = mappedRGB(:,2);
imgb(ind) = mappedRGB(:,3);

save('MappedMouseGeneVolumeSelectionSlice2D','imgr','imgg','imgb');

