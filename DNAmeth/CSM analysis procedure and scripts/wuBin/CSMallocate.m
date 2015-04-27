function g2 = CSMallocate(X, g1, delta)
% CSM detection step 1: (2) choose seed clusters and allocate reads to seed clusters
% X: matrix of size N by M, N: number of sites, M: number of reads
% g1: vector of length M, cluster ID, length(unique(g1)) > 0
% delta: scalar, threshold for choosing seed clusters, 0 < delta < 0.5
% g2: vector of length M, cluster ID, length(unique(g2)) = 2

[N, M] = size(X);
K = length(unique(g1));
c1 = zeros(1, K);
c2 = zeros(1, K);
c3 = zeros(1, K);
j1 = 0;
j2 = 0;
j3 = 0;
for i = 1 : K
    Xi = X(:, g1 == i);
	pih = mean(Xi, 2);
%     if(sum(pih < delta) == N)
    if(mean(pih) < delta)
        j1 = j1 + 1;
        c1(j1) = i;
    elseif(mean(pih) > 1 - delta)
        j2 = j2 + 1;
        c2(j2) = i;
    else
        j3 = j3 + 1;
        c3(j3) = i;
    end
end
if (sum(c1) == 0 || sum(c2) == 0)
    g2 = ones(1, M);
else
    g2 = zeros(1, M);
    L1 = c1(c1 ~= 0);
    L2 = c2(c2 ~= 0);
    L3 = c3(c3 ~= 0);
    ix_1 = ismember(g1, L1);
    g2(ix_1) = 1;
    X1 = X(:, ix_1);
    ix_2 = ismember(g1, L2);
    g2(ix_2) = 2;
    X2 = X(:, ix_2);
    ix_3 = ismember(g1, L3);
    X3 = X(:, ix_3);
    m3 = sum(ix_3);
    p1h_mat = repmat(mean(X1, 2), [1, m3]);
    p2h_mat = repmat(mean(X2, 2), [1, m3]);
    g2(ix_3) = (sum((X3 - p1h_mat) .^ 2) > sum((X3 - p2h_mat) .^ 2)) + 1;
end
