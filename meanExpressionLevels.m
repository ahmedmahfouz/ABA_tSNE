% Obtain mean expression levels by averaging probes

function meanExpression = meanExpressionLevels(dataDir, donor, probesPerGene)

% Determine which probes represent which gene
if nargin < 3
    probesPerGene = probeindices2genes(dataDir,donor);
end

% Load expression matrix
load([dataDir '/normalized_microarray_donor' donor '/MicroarrayExpression.mat']);

% For each gene
ngenes = length(probesPerGene);
nsamples = size(MicroarrayExpression,2);
meanExpression = zeros(ngenes,nsamples);
for gnr = 1:ngenes
    meanExpression(gnr,:) = mean(MicroarrayExpression(probesPerGene{gnr},:));
end
