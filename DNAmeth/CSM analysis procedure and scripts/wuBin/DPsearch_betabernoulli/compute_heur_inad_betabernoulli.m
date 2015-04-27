 
function l = compute_heur_inad_betabernoulli(s,marginals)
  N0 = length(s.c);
  l  = - marginals(N0+1); %zhu2012: note this is the sum of H(x_i) from N0+1 to N.
