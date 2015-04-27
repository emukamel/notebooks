% In this function we generate random samples from DP(alpha, G0)
% distribution following the Polye Urn scheme, where G0 is the baseline 
% distribution. Here we choose G0~ Dirichlet.

function [pX, c, K] = DP_rnd_beta(alpha, G0alphabeta, w, N)

% Input:
% alpha:         the precision parameter of DP distribution.
% G0alphabeta:   the hyperprior parameters in G0. G0alphabeta is 2 by w
%                matrix of alpha, beta parameters. 
% Reference: Page 45, "The Dirichlet process, related priors
% and posterior asymptotics", by Subhashis Ghosal.
% http://www4.stat.ncsu.edu/~ghoshal/papers/bayesasymp.pdf

% Output:
% X:  the data generated, N by w matrix.
% c:  the labels for cluster assignment.
% K:  the number of clusters in total.

pX  =  NaN(N, w);
c = NaN(N, 1);
pX(1,:) = betarnd(G0alphabeta(1,:), G0alphabeta(2,:)); 
c(1)   = 1;
K = 1;
for i = 2:N
    % Pr (Xi = Xj)    =  1/(alpha+i-1),j=1,...i-1
    % Pr (Xi = X_new) =  alpha/(alpha+i-1). 
    u = rand(1)*(alpha+i-1);
    csum = cumsum([1:(i-1),alpha]);
    l0 = find(csum - u>0);
    j = l0(1);
    if j <= (i-1) % this is the case Xi== Xj
       pX(i,:) = pX(j,:);        
       c(i)   = c(j);
    else % this is the case Xi = X_new
       pX(i,:) = betarnd(G0alphabeta(1,:), G0alphabeta(2,:));        
       K      = K + 1;
       c(i)   = K;
    end        
end
