% Jensen-Shannon divergence

function [d,histmatP,histmatQ] = jsdiv_pts(P,Q,nbins,minbin,maxbin)
    
[d,histmatP,histmatQ] = kldiv_pts(P,Q,nbins,minbin,maxbin);
d2 = kldiv_pts(Q,P,nbins,minbin,maxbin);

d = (d+d2)/2;
