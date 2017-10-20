%% Perform Barnes-Hut t-SNE on human gene data
clear variables
close all

addpath('../bh-tsne');
ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);

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
        selgenes = allgindices(indexdata.probe2gene_index>0);
    end
    
    % Load gene coordinates, color, structure_ID
    sampledata = load([datadir 'sample']);
    coords = [coords; [sampledata.sample.mni_x,sampledata.sample.mni_y,sampledata.sample.mni_z]];
    color_RGB = [color_RGB; sampledata.sample.color_RGB];
    structure_id = [structure_id; sampledata.sample.structure_id];
    
    brain_id = [brain_id; bnr*ones(length(sampledata.sample.mni_x),1)];
    
    % Normalize data
    Xnew = expressionlevels.MicroarrayExpression(selgenes,:)';
    meanXnew = mean(Xnew,2);
    stdXnew = std(Xnew,[],2);
    Xnew = bsxfun(@minus,Xnew,meanXnew);
    Xnew = bsxfun(@rdivide,Xnew,stdXnew);
    
    % Data subselection
    X = [X; Xnew];
end

%% Perform t-SNE
perplexity = 30;
theta = 0.7;
initial_dims = min(1000,size(X,2));
mappedX = fast_tsne(X,[],perplexity,theta);

%% Save mapped data
save(['MappedGenesGroupSelection' id],'X','mappedX',...
    'perplexity','theta',...
    'coords','color_RGB','structure_id','brain_id');

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

saveas(1,'MappedGenesGroupSelection');

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

saveas(2,'MappedGenesGroupSelectionRegions');

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

saveas(3,'MappedGenesGroupSelectionBrainID');

