%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');

% perplexity_vec = fliplr(2.^(9:10));
perplexity_vec = 256;
np = length(perplexity_vec);
initial_dims = 100;

%% Load data (currently for a single brain)
outdir = 'results_mouse_genes_perplexity_sweep_th0p2/';
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

%% Perform t-SNE
for pnr = 1:np
    
    X = expressionlevels.expressionMatrix;
    
    % selvec = randperm(size(X,1));
    sampleselvec = 1:size(X,1);
    geneselvec = selection.highConfGenes_coronal;
    X = X(sampleselvec,geneselvec);
    
    % Normalize
    X = zscore(X);
    
    coords = coords(sampleselvec,:);
    theta = 0.2;
	perplexity = perplexity_vec(pnr);
    
    fprintf('Perplexity %d (%d/%d)\n',perplexity,pnr,np);
    mappedX = fast_tsne(X,2,initial_dims,perplexity,theta);
    
    %
    minx = min(mappedX);
    maxx = max(mappedX);
    maxval = max([abs(minx) abs(maxx)]);
    mappedab = mappedX*127/maxval;
    tC = makecform('lab2srgb');
    mappedRGB = applycform([50*ones(size(mappedab,1),1) mappedab],tC);
    
    %% Save mapped data
    save([outdir sprintf('MappedMouseGenesPerplexitySweep2Dperp%did100',perplexity)],'X','mappedX','coords',...
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
    
    saveas(3,[outdir sprintf('MappedMousePerplexitySweep2Dperp%did100',perplexity)]);
    saveas(3,[outdir sprintf('MappedMousePerplexitySweep2Dperp%did100.png',perplexity)]);
    
    %% Image volume
    isz = [67 41 58];
    imgr = zeros(isz);
    imgg = zeros(isz);
    imgb = zeros(isz);
    
    ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
    imgr(ind) = mappedRGB(:,1);
    imgg(ind) = mappedRGB(:,2);
    imgb(ind) = mappedRGB(:,3);
    
    save([outdir sprintf('MappedMouseGeneVolumePerplexitySweep2Dperp%did100',perplexity)],'imgr','imgg','imgb');
end
