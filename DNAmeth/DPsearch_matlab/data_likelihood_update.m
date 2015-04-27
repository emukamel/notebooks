function l = data_likelihood_update(X,alpha,c,phi,cnt)
  [N,W] = size(X);
  M     = length(c);
  
  if M <= 1,
    l = data_likelihood(X,alpha,c);
  else
    k = c(M);
    
    if cnt(k) == 1, % new cluster
      l = data_likelihood(X(M,:),alpha,1);
    else            % existing cluster
      xm = X(M,:);
      sx = phi(k,:);
      
      l = gammaln(sum(xm)+1) - sum(gammaln(xm+1)) + ...
          sum(gammaln(alpha+sx)) - gammaln(sum(alpha+sx)) + ...
          gammaln(sum(alpha+sx-xm)) - sum(gammaln(alpha+sx-xm));
    end;
  end;
