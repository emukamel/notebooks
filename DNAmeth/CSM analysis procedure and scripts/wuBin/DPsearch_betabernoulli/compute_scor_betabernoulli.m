function [l,dl] = compute_scor_betabernoulli(s,X,alpha,G0alphabeta)
  N0 = length(s.c);
  if N0 == 1
    dl = data_likelihood_betabernoulli(X,G0alphabeta,s.c);
  else
    dl = s.par.l + data_likelihood_update_betabernoulli(X,G0alphabeta,s.c,s.phi,s.cnt);
  end 
  l  = - dl - log_DP_prior_count_complete2(s.c, alpha, size(X,1), s.m); 
  % zhu: negative log likelihood H(x_1)p(c)
