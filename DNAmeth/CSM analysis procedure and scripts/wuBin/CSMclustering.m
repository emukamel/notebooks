function g = CSMclustering(Y)
% CSM detection step 1: (1) clustering data (need to include several different clustering methods)
% Y: matrix of size M by N, M: number of reads, N: number of sites
% g: vector of length M, predicted cluster ID

[M, N] = size(Y);

alpha = 1; 
G0alphabeta = 0.1 .* ones(2, N);
beamSize = 20; % 20-100
[r, ~] = DPsearch_betabernoulli(Y, alpha, G0alphabeta, beamSize, ones(1, M), 'i');
g = r.c;
