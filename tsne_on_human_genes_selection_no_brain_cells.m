%% Perform Barnes-Hut t-SNE on human gene data
clear variables
close all

addpath('../bh-tsne');
ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);
celltypes = {'neurons','oligodendrocytes','astrocytes'};
ncelltypes = length(celltypes);

for bnr = 1:nbrains
%% Load data (currently for a single brain)
id = ids{bnr};
datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
expressionlevels = load([datadir 'MicroarrayExpression']);

% Load gene selection
indexdata = load([datadir 'probe2gene_index']);
celldata = load([datadir 'cell_type_genes_indices_old']);
allgindices = 1:size(expressionlevels.MicroarrayExpression,1);
selnobraincells = true(size(expressionlevels.MicroarrayExpression,1),1);
for cnr = 1:ncelltypes
    selnobraincells(celldata.geneSubSets.(celltypes{cnr})) = false;
end
selgenes = allgindices(indexdata.probe2gene_index>0 & selnobraincells);
ngenes = length(selgenes);

% Load gene coordinates
coordsdata = load([datadir 'sample']);
coords = [coordsdata.sample.mni_x,coordsdata.sample.mni_y,coordsdata.sample.mni_z];

%% Perform t-SNE
% Data subselection
X = expressionlevels.MicroarrayExpression';
X = X(:,selgenes);

% Perform t-SNE
perplexity = 30;
theta = 0.7;
initial_dims = min(1000,size(X,2));
mappedX = fast_tsne(X,initial_dims,perplexity,theta);

%% Save mapped data
save(['MappedGenesSelectionNoBrainCells' id],'X','mappedX','selgenes','coords',...
    'initial_dims','perplexity','theta');

%% Visualize mapped data
% Colorspace
crange = max(coords)-min(coords);
mincoords = min(coords);
nsamples = size(X,1);

%
figure(1)
clf
hold on
for snr = 1:nsamples
    c = (coords(snr,:)-mincoords)./crange;
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
end

% saveas(1,'MappedGenesSelection');
end