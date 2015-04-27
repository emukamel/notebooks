function [pval, d] = CSMtest_chisq(X, g, tau)
% CSM detection step 2: test hypothesis for clustered data (Wald test)
% X: matrix of size N by M, N: number of sites, M: number of reads
% g: vector of length M, cluster ID
% tau: scalar, separation parameter
% pval: testing p-value

X1 = X(:, g == 1);
m1 = size(X1, 2);
p1h = mean(X1, 2);
X2 = X(:, g == 2);
m2 = size(X2, 2);
p2h = mean(X2, 2);
ph = mean(X, 2);
if (sum(ph == 0) == 0 && sum(ph == 1) == 0)
    if (sum(p1h >= p2h) == 4) || (sum(p1h <= p2h) == 4)
        s = min((abs(p1h - p2h) - tau) .^ 2 ./ (p1h .* (1 - p1h) ./ m1 + p2h .* (1 - p2h) ./ m2));
        pval = 1 - chi2cdf(s, 1);
    else
        pval = 1.1;
    end
else
	pval = []; % NOT POLYMORPHIC!
end
d = abs(mean(p1h) - mean(p2h));
