%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');

% initial_dims_vec = [2 3 5 10:10:100];
% initial_dims_vec = [2 3 5 10 15 20];
initial_dims_vec = 10;
nids = length(initial_dims_vec);

%% Load data (currently for a single brain)
outdir = 'results_mouse_genes_selection_cell_types/';
datadir = '/home/mvandegiessen/data/tSNE_ABA/AllenMouseBrain_sagittal/';
expressionlevels = load([datadir 'expressionMatrix']);
ngenes = size(expressionlevels.expressionMatrix,2);

%% Load gene selection
selection = load([datadir 'highConfGenesInd_sagittal']);
selection.highConfGenes_sagittal = selection.highConfGenes_sagittal-4345;
cellselection = load([datadir 'cell_type_genes_indices_sagittal']);
colors = load([datadir 'voxels']);

% Load gene coordinates
coordsdata = load([datadir 'voxels']);
[cx,cy,cz] = ind2sub([67 41 58],coordsdata.voxel.brain_idx);
coords = [cx cy cz];

%% Perform t-SNE
for inr = 1:nids
    X = expressionlevels.expressionMatrix;
    
    % selvec = randperm(size(X,1));
    sampleselvec = 1:size(X,1);
    geneselvec = selection.highConfGenes_sagittal(cellselection.geneSubSets.oligodendrocytes);
    X = X(sampleselvec,geneselvec);
    
    % Normalize
    X = zscore(X);
    
    coords = coords(sampleselvec,:);
    perplexity = 30;
    theta = 0.7;
    initial_dims = min(initial_dims_vec(inr),size(X,2));
    mappedX = fast_tsne(X,2,initial_dims,perplexity,theta);
    
    % Mapped colors
    minx = min(mappedX);
    maxx = max(mappedX);
    maxval = max([abs(minx) abs(maxx)]);
    mappedab = mappedX*127/maxval;
    tC = makecform('lab2srgb');
    mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);
    
    %% Save mapped data
    save([outdir sprintf('MappedMouseGenesSelection2DsagittalOligodendrocytesid%d',initial_dims)],'X','mappedX','coords',...
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
        if ~isempty(colors.voxel.color_HEX{snr})
            colorRGB = rgbconv(colors.voxel.color_HEX{snr});
        else
            colorRGB = [0 0 0];
        end
        plot(mappedX(snr,1),mappedX(snr,2),'.','Color',colorRGB);
    end
    
    saveas(3,[outdir sprintf('MappedMouseGenesSelection2DsagittalOligodendrocytes%d',initial_dims)]);
    saveas(3,[outdir sprintf('MappedMouseGenesSelection2DsagittalOligodendrocytes%d.png',initial_dims)]);
    
    %% Image volume
    isz = [67 41 58];
    imgr = zeros(isz);
    imgg = zeros(isz);
    imgb = zeros(isz);
    
    ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
    imgr(ind) = mappedRGB(:,1);
    imgg(ind) = mappedRGB(:,2);
    imgb(ind) = mappedRGB(:,3);
    
    save([outdir sprintf('MappedMouseGeneVolumeSelection2DsagittalOligodendrocytesid%d',initial_dims)],'imgr','imgg','imgb');
end
