%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');

%% Load data (currently for a single brain)
% id = '9861';
% id = '10021';
id = 'all';
ids = {'9861','10021'};

if strcmp(id,'all')
    X = [];
    ontologycolor = [];
    for inr = 1:length(ids)
        datadir = ['/home/mvandegiessen/data/tSNE_ABA/rnaSeq_21April2014/rnaseq_donor' ids{inr} '/'];
        expressionlevels = load([datadir 'rnaSequencing']);
        X = [X; expressionlevels.rnaSequencing'];
        
        load([datadir 'sample']);
        ontologycolor = [ontologycolor; sample.ontology_color];
    end
    ngenes = size(expressionlevels.rnaSequencing,1);
else
    datadir = ['/home/mvandegiessen/data/tSNE_ABA/rnaSeq_21April2014/rnaseq_donor' id '/'];
    expressionlevels = load([datadir 'rnaSequencing']);
    X = expressionlevels.rnaSequencing';
    ngenes = size(expressionlevels.rnaSequencing,1);
    
    load([datadir 'sample']);
    ontologycolor = sample.ontology_color;
end

%% Perform t-SNE

% selvec = randperm(size(X,1));
nsamples = size(X,1);
selvec = 1:nsamples;
X = X(selvec,:);
perplexity = 30;
theta = 0.1;
initial_dims = min(300,size(X,2));
mappedX = fast_tsne(X,2,initial_dims,perplexity,theta);

%% Save mapped data
save(['MappedRNAseq2D' id],'X','mappedX',...
    'initial_dims','perplexity','theta','selvec');

%% Visualize mapped data based on t-SNE coordinates
minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);
mappedab = mappedX*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);

figure(3)
clf
hold on
for snr = 1:nsamples
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB(snr,:));
end

saveas(3,['MappedRNAseqLabMap2D' id]);
saveas(3,['MappedRNAseqLabMap2D' id '.png']);

%% Color based visualization

figure(4)
clf
hold on
for snr = 1:nsamples
    tmp = ontologycolor{snr};
    c = rgbconv(tmp(2:7));
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
end

saveas(4,['MappedRNAseqRegion2D' id]);
saveas(4,['MappedRNAseqRegion2D' id '.png']);
