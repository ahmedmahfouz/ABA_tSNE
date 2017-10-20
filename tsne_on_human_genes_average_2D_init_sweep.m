%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bhtsne_intermediate');
ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);

% initial_dims_vec = [2 3 5 10:10:100];
initial_dims_vec = [2 3 4 5 6];
% initial_dims_vec = 300;
nids = length(initial_dims_vec);

%% Load data
for bnr = 1:1%nbrains
    id = ids{bnr};
    datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
    outdir = 'results_human_average_pca_sweep_fixed_init/';
    expressionlevels = load([datadir 'geneMicroarrayExpression']);
    
    %     % Load gene selection
    %     if bnr == 1
    %         indexdata = load([datadir 'probe2gene_index']);
    %         allgindices = 1:size(expressionlevels.MicroarrayExpression,1);
    %         selprobes = allgindices(indexdata.probe2gene_index>0);
    %     end
    
    % Load gene coordinates, color, structure_ID
    sampledata = load([datadir 'sample']);
    coords = [sampledata.sample.mni_x,sampledata.sample.mni_y,sampledata.sample.mni_z];
    color_RGB = sampledata.sample.color_RGB;
    structure_id = sampledata.sample.structure_id;
    brain_id = bnr*ones(length(sampledata.sample.mni_x),1);
    
    % Normalize data
    X0 = expressionlevels.geneMicroarrayExpression';
    X0 = zscore(X0,0,1);
    
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
        no_dims = 2;
        randomseed = 0;
%         mappedX = fast_tsne(X,2,initial_dims,perplexity,theta);
        
        save_interval = 50;
        if inr == 1
            mappedX = fast_tsne_intermediate(X, no_dims, initial_dims, perplexity, theta, randomseed, [], save_interval);
        else
            mappedX = fast_tsne_intermediate(X, no_dims, initial_dims, perplexity, theta, randomseed, mappedX, save_interval);
        end
        
        destdir = sprintf('initial_dims=%04d',initial_dims);
        if ~exist([outdir destdir],'dir')
            mkdir([outdir destdir]);
        end
        movefile('intermediate*.dat',[outdir destdir]);
%         files = dir('intermediate*.dat');
%         for i=1:length(files)
%             mappedX = read_data(files(i).name);
%             plot(mappedX(:,1), mappedX(:,2),'.');
%             pause(1);
%         end
        
        %% Save mapped data
        save([outdir destdir '/' sprintf('MappedHumanGenesAverage2D%sid%d',id,initial_dims)],'X','mappedX','coords',...
            'initial_dims','perplexity','theta','selvec','color_RGB','structure_id','brain_id');
        
        %% Visualize mapped data based on t-SNE coordinates
        minx = min(mappedX);
        maxx = max(mappedX);
        maxval = max([abs(minx) abs(maxx)]);
        mappedab = mappedX*127/maxval;
        tC = makecform('lab2srgb');
        mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);
        
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
        
        saveas(2,[outdir destdir '/' sprintf('MappedHumanGenesAverageRegions2D%sid%d',id,initial_dims)]);
    end
end
