%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');
ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);

% initial_dims_vec = [2 3 5 10:10:100];
initial_dims_vec = [1 4 6];
% initial_dims_vec = 11:19;
nids = length(initial_dims_vec);

%% Load data
X = [];
coords = [];
color_RGB = [];
structure_id = [];
brain_id = [];
for bnr = 1:nbrains
    id = ids{bnr};
    datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
    expressionlevels = load([datadir 'MicroarrayExpression']);
    
    % Load gene selection
    if bnr == 1
        indexdata = load([datadir 'probe2gene_index']);
        allgindices = 1:size(expressionlevels.MicroarrayExpression,1);
        selprobes = allgindices(indexdata.probe2gene_index>0);
    end
    
    % Load gene coordinates, color, structure_ID
    sampledata = load([datadir 'sample']);
    coords = [coords; [sampledata.sample.mni_x,sampledata.sample.mni_y,sampledata.sample.mni_z]];
    color_RGB = [color_RGB; sampledata.sample.color_RGB];
    structure_id = [structure_id; sampledata.sample.structure_id];
    
    brain_id = [brain_id; bnr*ones(length(sampledata.sample.mni_x),1)];
    
    % Normalize data
    Xnew = expressionlevels.MicroarrayExpression(selprobes,:)';
%     meanXnew = mean(Xnew,2);
%     stdXnew = std(Xnew,[],2);
%     Xnew = bsxfun(@minus,Xnew,meanXnew);
%     Xnew = bsxfun(@rdivide,Xnew,stdXnew);
    
    % Data subselection
    X = [X; Xnew];
end
X0 = X;

%% Perform t-SNE
for inr = 1:nids
    X = X0;
    
    % selvec = randperm(size(X,1));
    selvec = 1:size(X,1);
    X = X(selvec,:);
    coords = coords(selvec,:);
    perplexity = 30;
    theta = 0.5;
    initial_dims = min(initial_dims_vec(inr),size(X,2));
    mappedX = fast_tsne(X,2,initial_dims,perplexity,theta);
    
    %% Save mapped data
    save(sprintf('MappedHumanGenesSelection2Did%d',initial_dims),'X','mappedX','coords',...
        'initial_dims','perplexity','theta','selvec','selprobes','color_RGB','structure_id','brain_id');
    
    
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
    minx = min(mappedX);
    maxx = max(mappedX);
    maxval = max([abs(minx) abs(maxx)]);
    mappedab = mappedX*127/maxval;
    tC = makecform('lab2srgb');
    mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);
    
    % figure(3)
    % clf
    % hold on
    % for snr = 1:nsamples
    %     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB(snr,:));
    % end
    %
    % saveas(3,'MappedMouseGenesSelectionLabMap2D300');
    % saveas(3,'MappedMouseGenesSelectionLabMap2D300.png');
    
    %% Image volume
%     isz = [67 41 58];
%     imgr = zeros(isz);
%     imgg = zeros(isz);
%     imgb = zeros(isz);
%     
%     ind = sub2ind(isz,coords(:,1),coords(:,2),coords(:,3));
%     imgr(ind) = mappedRGB(:,1);
%     imgg(ind) = mappedRGB(:,2);
%     imgb(ind) = mappedRGB(:,3);
%     
%     save(sprintf('MappedHumanGeneVolumeSelection2Did%d',initial_dims),'imgr','imgg','imgb');
    %% Visualize mapped data using region colors
    nsamples = size(X,1);
    
    %
    figure(2)
    clf
    hold on
    for snr = 1:nsamples
        c = color_RGB(snr,:);
        plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
    end
    
    saveas(2,sprintf('MappedHumanGenesGroupSelectionRegions2Did%d',initial_dims));
    
    %% Visualize mapped data using brain id colors
    nsamples = size(X,1);
    
    %
    figure(3)
    clf
    hold on
    cm = hsv(nbrains);
    for bnr = 1:nbrains
        c = cm(bnr,:);
        plot(mappedX(brain_id==bnr,1),mappedX(brain_id==bnr,2),'.','Color',c);
    end
    
    saveas(3,sprintf('MappedHumanGenesGroupSelectionBrainID2Did%d',initial_dims));
end
