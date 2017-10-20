%%% 25 Feb 2014
%%% collapse probes to genes 

function probe2gene_index = probes2genes(dataDir, donor)

% load the donor probe info and the expression matrix
load([dataDir '/normalized_microarray_donor' donor '/probe.mat']);
load([dataDir '/normalized_microarray_donor' donor '/MicroarrayExpression.mat']);

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
probe2gene_index = zeros(length(probe.probe_id),1,'int8');
geneVar = var(MicroarrayExpression');
uniqueGenes = unique(probe.gene_symbol(intersect(probes_entrezId,probes_symbol)));
for g = 1 : length(uniqueGenes)
    currGeneInd = find(strcmpi(probe.gene_symbol,uniqueGenes{g})==1);
    if length(currGeneInd) > 1
        [maxVar maxVarInd] = max(geneVar(currGeneInd));
        probe2gene_index(currGeneInd(maxVarInd)) = 1;
    else
        probe2gene_index(currGeneInd) = 1;
    end
end










