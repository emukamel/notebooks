function m = counts(c,N)  % this assumes that there are no "empty"
                          % clusters ... safe for search, unsafe for
                          % sampling (unless you GC every iteration)
  k = max(c);
  m = zeros(1,N);
  t = zeros(1,k);
  for i=1:length(c), t(c(i)) = t(c(i)) + 1; end;
  for i=1:k, m(t(i)) = m(t(i)) + 1; end;
  