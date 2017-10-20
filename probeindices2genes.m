% Select probe indices that belong to a gene


% Based on probes2genes
%%% 25 Feb 2014
%%% collapse probes to genes 

function probe2gene_indices = probeindices2genes(dataDir, donor)

% load the donor probe info and the expression matrix
load([dataDir '/normalized_microarray_donor' donor '/probe.mat']);
% load([dataDir '/normalized_microarray_donor' donor '/MicroarrayExpression.mat']);

%% exclude probes with no entrez id
probes_entrezId = find(probe.entrez_id ~= 0);
%% exclude probes with no symbol
probes_symbol =  find(strcmpi(probe.gene_symbol, 'na') ~= 1);

%% for each gene, select the "best" probe (LATER)
% NOTE:  REQUIRES CALCULATING THE FULL CORRELATION MATRIX
% silmilar to the collapseRows function in the WCGNA package. (1) for genes 
% with 3 or more probes: select the probe with the max connectivity and (2)
% for genes with 2 probes: select the probe with the higher variance

%% for each gene, select the "best" probe
% (ALTERNATIVE) select the probe with the highest variance
uniqueGenes = unique(probe.gene_symbol(intersect(probes_entrezId,probes_symbol)));
probe2gene_indices = cell(length(uniqueGenes),1);
for g = 1 : length(uniqueGenes)
    currGeneInd = find(strcmpi(probe.gene_symbol,uniqueGenes{g})==1);
    probe2gene_indices{g} = currGeneInd;
end










