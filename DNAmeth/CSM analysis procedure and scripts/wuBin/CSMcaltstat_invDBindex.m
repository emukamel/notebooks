function tstat = CSMcaltstat_invDBindex(X, g)
% calculate test statistic for permutation test by Davies-Bouldin index
% X: matrix of size N by M, N: number of sites, M: number of reads
% g: vector of length M, cluster ID

N = size(X, 1);
K = length(unique(g));
pih = NaN(N, K);
intradist = zeros(1, K);
for i = 1 : K
	Xi = X(:, g == i);
    mi = sum(g == i);
	centroidi = mean(Xi, 2);
    pih(:, i) = centroidi;
    intradist(i) = mean(sqrt(sum((Xi - repmat(centroidi, [1, mi])) .^ 2)));
end
pairmat = combnk(1 : K, 2);
npairs = size(pairmat, 1);
interdist = zeros(1, npairs);
pairdist = zeros(1, npairs);
for i = 1 : npairs
    cid1 = pairmat(i, 1);
    cid2 = pairmat(i, 2);
    interdist(i) = sqrt(sum((pih(:, cid1) - pih(:, cid2)) .^ 2));
    pairdist(i) = interdist(i) / (intradist(cid1) + intradist(cid2));
end
tstat = sum(pairdist);
