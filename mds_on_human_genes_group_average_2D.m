%% Perform Barnes-Hut t-SNE on mouse gene data
clear variables
close all

addpath('../bh-tsne');
ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);
outdir = 'results_human_average_group_MDS/';

%% Load data
X = [];
coords = [];
color_RGB = [];
structure_id = [];
brain_id = [];
for bnr = 1:nbrains
    id = ids{bnr};
    datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
%     expressionlevels = load([datadir 'MicroarrayExpression']);
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
    % Per gene
    Xnew = zscore(Xnew,0,1);
    
    % Data subselection
    X = [X; Xnew];
end
X0 = X;

%% Compute distance matrix
disp('Distance matrix');
% D = pdist(X,'euclidean');
nsamples = size(X,1);
D = zeros(nsamples,nsamples);
tic
parfor snr1 = 1:nsamples
%     fprintf('Sample %d/%d\n',snr1,nsamples);
    dist = zeros(1,nsamples);
    for snr2 = 1:nsamples
        dist(snr2) = sqrt(1-xcov(X(snr1,:),X(snr2,:),0,'coeff')^2);
    end
    D(snr1,:) = dist;
end
toc

D = abs(D);
D(eye(size(D))>0) = 0;

%% MDS
disp('2D MDS');
mappedX2 = mdscale(D,2);
disp('3D MDS');
mappedX3 = mdscale(D,3);

%% Save mapped data
save([outdir 'MDSHumanGenesGroupAverageMeanSub'],'D','X','mappedX2','mappedX3','coords','color_RGB','structure_id','brain_id');

%% Visualize mapped data based on SVD coordinates
nsamples = size(X,1);

figure(2)
clf
hold on
for snr = 1:nsamples
    c = color_RGB(snr,:);
    plot3(mappedX3(snr,1),mappedX3(snr,2),mappedX3(snr,3),'.','Color',c);
end
view(3)
axis equal

saveas(2,[outdir 'MDSHumanGenesGroupAverageMeanSubRegions3D']);

%% Visualize mapped data using brain id colors
nsamples = size(X,1);

%
figure(3)
clf
hold on
cm = hsv(nbrains);
for bnr = 1:nbrains
    c = cm(bnr,:);
    plot3(mappedX3(brain_id==bnr,1),mappedX3(brain_id==bnr,2),mappedX3(brain_id==bnr,3),'.','Color',c);
end
view(3)
axis equal

saveas(3,[outdir 'MDSHumanGenesGroupAverageMeanSubBrainID3D']);

%% Visualize mapped data based on SVD coordinates
figure(4)
clf
hold on
for snr = 1:nsamples
    c = color_RGB(snr,:);
    plot(mappedX2(snr,1),mappedX2(snr,2),'.','Color',c);
end

saveas(4,[outdir 'MDSHumanGenesGroupAverageMeanSubRegions2D']);

%% Visualize mapped data using brain id colors
nsamples = size(X,1);

%
figure(5)
clf
hold on
cm = hsv(nbrains);
for bnr = 1:nbrains
    c = cm(bnr,:);
    plot(mappedX2(brain_id==bnr,1),mappedX2(brain_id==bnr,2),'.','Color',c);
end

saveas(5,[outdir 'MDSHumanGenesGroupAverageMeanSubBrainID2D']);
