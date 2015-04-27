function [l, marg_like]= log_likelihood_given_c(data, c, G0alpha)
% This function is written by zhu to test whether the log likelihood
% calcualted using the DPsearch algorithm is equal the true. 
% This function calculate the loglikelihood give a known clustering.
% data: X.
% c: a vector of cluster asignment. 

% Output:
% marg_like: the vector of k, each is the marginal likelihood of each
% cluster.
% l: the sum of marg_like

k = length(unique(c));
marg_like = NaN(k, 1);

for i = 1:k
     marg_like(i)= data_likelihood_betabernoulli(data(c==i,:),G0alpha,c(c==i));
end

l = sum(marg_like);