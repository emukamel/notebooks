function pval = CSMtest_perm(X, g, nperm)
% CSM detection step 2: test hypothesis for clustered data (permutation test)
% X: matrix of size N by M, N: number of sites, M: number of reads
% g: vector of length M, cluster ID
% nperm: number of random permutations to get null distribution
% pval: testing p-value

M = size(X, 2);
% tstatobs = CSMcaltstat_interdist(X, g);
tstatobs = CSMcaltstat_invDBindex(X, g);
tstatdist = NaN(1, nperm);
for i = 1 : nperm
    gperm = randsample(g, M);
%     tstatdist(i) = CSMcaltstat_interdist(X, gperm);
    tstatdist(i) = CSMcaltstat_invDBindex(X, gperm);
end
pval = sum(tstatdist > tstatobs) / nperm;
