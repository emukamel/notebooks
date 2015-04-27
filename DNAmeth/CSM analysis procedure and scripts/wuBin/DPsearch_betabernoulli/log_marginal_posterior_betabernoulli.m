function l = log_marginal_posterior_betabernoulli(x, G0alphabeta, xs)
  N = size(xs,1);

  if N > 0,
    l = data_likelihood_betabernoulli([x;xs], G0alphabeta, ones(1,N+1)) - ...
        data_likelihood_betabernoulli(   xs , G0alphabeta, ones(1,N  ));
  else
    l = data_likelihood_betabernoulli(x, G0alphabeta, 1);
  end;
  