%% Perform Barnes-Hut t-SNE on human gene data
clear variables
close all

addpath('../bh-tsne');

%% Load data (currently for a single brain)
id = '15697';
datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
expressionlevels = load([datadir 'MicroarrayExpression']);
ngenes = size(expressionlevels.MicroarrayExpression,1);

% Load gene coordinates
coordsdata = load([datadir 'sample']);
coords = [coordsdata.sample.mni_x,coordsdata.sample.mni_y,coordsdata.sample.mni_z];

%% Perform t-SNE
X = expressionlevels.MicroarrayExpression';
perplexity = 30;
theta = 0.7;
initial_dims = min(1000,size(X,2));
mappedX = fast_tsne(X,initial_dims,perplexity,theta);

%% Save mapped data
save(['MappedGenes' id],'X','mappedX','coords',...
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

saveas(1,['MappedGenes' id]);
