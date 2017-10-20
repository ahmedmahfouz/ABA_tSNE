%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');
ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);

% initial_dims_vec = [2 3 5 10:10:100];
initial_dims_vec = [1 2 3 4 5 6 10:10:100];
initial_dims_vec = [7 8 9];
% initial_dims_vec = 100;
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
    outdir = 'results_human_average_pca_sweep/';
    expressionlevels = load([datadir 'geneMicroarrayExpression']);
    
    %     % Load gene selection
    %     if bnr == 1
    %         indexdata = load([datadir 'probe2gene_index']);
    %         allgindices = 1:size(expressionlevels.MicroarrayExpression,1);
    %         selprobes = allgindices(indexdata.probe2gene_index>0);
    %     end
    
    % Load gene coordinates, color, structure_ID
    sampledata = load([datadir 'sample']);
    coords = [coords; [sampledata.sample.mni_x,sampledata.sample.mni_y,sampledata.sample.mni_z]];
    color_RGB = [color_RGB; sampledata.sample.color_RGB];
    structure_id = [structure_id; sampledata.sample.structure_id];
    brain_id = [brain_id; bnr*ones(length(sampledata.sample.mni_x),1)];
    
    % Normalize data
    Xnew = expressionlevels.geneMicroarrayExpression';
    %     meanXnew = mean(Xnew,2);
    %     stdXnew = std(Xnew,[],2);
    %     Xnew = bsxfun(@minus,Xnew,meanXnew);
    %     Xnew = bsxfun(@rdivide,Xnew,stdXnew);
    % Per gene
%     Xnew = zscore(Xnew,0,1);
    % Per brain
%     Xnew = zscore(Xnew,0,2);
    
    % Data subselection
    X = [X; Xnew];
end
X0 = X;

%% Compute mean of absolute differences
madX = mad(X,1,1);
mmadX = madX./median(X);

%% Show histogram
[n,x] = hist(madX,500);
figure(1)
clf
bar(x,n);

[mn,mx] = hist(mmadX,500);
figure(2)
clf
bar(mx,mn);

%% Select data
thmadX = 0.154;
selvec = madX>thmadX;

[mn2,mx2] = hist(madX(selvec)./median(X(:,selvec)),500);
figure(3)
clf
bar(mx2,mn2);

save([datadir '../' 'mad_human_gene_selection'],'selvec','thmadX');

