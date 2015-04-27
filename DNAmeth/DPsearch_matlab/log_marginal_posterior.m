function l = log_marginal_posterior(x, alpha, xs)
  N = size(xs,1);

  if N > 0,
    l = data_likelihood([x;xs], alpha, ones(1,N+1)) - ...
        data_likelihood(   xs , alpha, ones(1,N  ));
  else
    l = data_likelihood(x, alpha, 1);
  end;
  