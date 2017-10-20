% Compute KL divergence for point clouds
%
% Should be extended to N-dimensional, currently 2-dimensional

function [d,histmatP,histmatQ] = kldiv_pts(P,Q,nbins,minbin,maxbin)

% Number of dimensions
ndims = size(P,2);

% Histogram allocations
if length(nbins) == 1
    nbins = repmat(nbins,[1 ndims]);
end
histmatP = zeros(nbins);
histmatQ = zeros(nbins);

% Bin boundaries
if nargin < 4
    minbin = prctile([P;Q],1);
end
if nargin < 5
    maxbin = prctile([P;Q],99);
end
binedges = cell(1,ndims);
binwidth = zeros(1,ndims);
for dnr = 1:ndims
    binwidth(dnr) = (maxbin(dnr)-minbin(dnr))/nbins(dnr);
    binedges{dnr} = minbin(dnr):binwidth(dnr):maxbin(dnr);
end

% Put all values within the boundaries
for dnr = 1:ndims
    selvec = P(:,dnr)<minbin(dnr);
    P(selvec,dnr) = minbin(dnr);
    selvec = P(:,dnr)>=maxbin(dnr)-binwidth(dnr)/2;
    P(selvec,dnr) = maxbin(dnr)-binwidth(dnr)/2;
    
    selvec = Q(:,dnr)<minbin(dnr);
    Q(selvec,dnr) = minbin(dnr);
    selvec = Q(:,dnr)>=maxbin(dnr)-binwidth(dnr)/2;
    Q(selvec,dnr) = maxbin(dnr)-binwidth(dnr)/2;
end

% Histogram
sP = floor((P-repmat(minbin,[size(P,1) 1]))./repmat(binwidth,[size(P,1) 1])+1);
sQ = floor((Q-repmat(minbin,[size(Q,1) 1]))./repmat(binwidth,[size(Q,1) 1])+1);

sPind = sub2ind(size(histmatP),sP(:,1),sP(:,2));
sQind = sub2ind(size(histmatQ),sQ(:,1),sQ(:,2));

hP = hist(sPind,1:numel(histmatP));
hQ = hist(sQind,1:numel(histmatQ));

histmatP = reshape(hP,size(histmatP));
histmatQ = reshape(hQ,size(histmatQ));

% Divergence
d = kldiv(eps+hP/sum(hP(:)),eps+hQ/sum(hQ(:)));
