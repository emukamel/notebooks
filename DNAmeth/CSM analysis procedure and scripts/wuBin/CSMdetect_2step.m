function [pval, gh2, d] = CSMdetect_2step(X, g, delta, tau, nperm)
% detect CSM using 2-step method (DPsearch and permutation test)
% X: matrix of size N by M, N: number of sites, M: number of reads
% g: row vector of length M, true cluster ID, optional
% pval: testing p-value

% STEP 1: Clustering
gh1 = CSMclustering(X');
gh2 = CSMallocate(X, gh1, delta);
d = 0;

% STEP 2: Testing
if length(unique(gh2)) == 1
% if length(unique(g)) == 1 % to evaluate the performance of permutation test only
    pval = 1;
else
%    pval = CSMtest_perm(X, gh2, nperm);
%     pval = CSMtest_perm(X, g, nperm); % to evaluate the performance of permutation test only
	[pval, d] = CSMtest_chisq(X, gh2, tau);
end

