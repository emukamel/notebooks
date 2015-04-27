function l = data_likelihood_update_betabernoulli(X,G0alphabeta,c,phi,cnt)
  %[,W] = size(X);
  M     = length(c);
  
  if M <= 1,
    l = data_likelihood_betabernoulli(X,G0alphabeta,c);
  else
    k = c(M);
    
    if cnt(k) == 1, % new cluster
      l = data_likelihood_betabernoulli(X(M,:),G0alphabeta,1);
    else            % existing cluster
      xm = X(M,:);
      sx = phi(k,:);
      
%       l = gammaln(sum(xm)+1) - sum(gammaln(xm+1)) + ...
%           sum(gammaln(alpha+sx)) - gammaln(sum(alpha+sx)) + ...
%           gammaln(sum(alpha+sx-xm)) - sum(gammaln(alpha+sx-xm));
     l = sum(gammaln(sx + G0alphabeta(1,:)) + gammaln(cnt(k)-sx+G0alphabeta(2,:))-...
         gammaln(sx-xm + G0alphabeta(1,:))- gammaln(cnt(k)-1 -(sx-xm)+G0alphabeta(2,:))+...
         gammaln(sum(G0alphabeta,1)+cnt(k)-1)-gammaln(sum(G0alphabeta,1)+cnt(k)));
    end;
  end;
