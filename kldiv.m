% 2D Kullback-Leibler divergence

function d = kldiv2d(P,Q)

% Compute entropy
e = log(P./Q).*P;

% Sum
d = sum(e(:));
