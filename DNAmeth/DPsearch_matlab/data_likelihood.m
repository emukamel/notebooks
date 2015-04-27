function l = data_likelihood(X,alpha,c)
  [N,W] = size(X);
  if length(c) > 1,
    fail;
  end;
  
  l = gammaln(W*alpha) - W*gammaln(alpha) + ...
      sum(gammaln(alpha + X(1,:))) - gammaln(sum(alpha + X(1,:))) + ...
      gammaln(sum(X(1,:)) + 1) - sum(gammaln(X(1,:)+1));
  
    
  
%   l = 0;
%   for i=1:length(K),
%     k  = K(i);
%     ck = find(c==k);
    
%     l = l + gammaln(W*alpha) - W * gammaln(alpha);
%     p = zeros(1,W);
%     for i=1:length(ck),
%       x = X(ck(i),:);

%       l = l + gammaln(sum(x) + 1) - sum(gammaln(x + 1));
%       p = p + x;
%     end;
%     l = l + sum(gammaln(alpha + p)) - gammaln(W*alpha + sum(p));
%   end;
  
  
%  a = alpha * W;
%  
%  l = 0;
%  for i=1:length(K),
%    k  = K(i);
%    ck = find(c==k);
%
%    l = l + gammaln(a) - W * gammaln(alpha);
%    p = ones(1,W) * alpha;
%    for i=1:length(ck),
%      x = X(ck(i),:);
%      p = p + x;
%      l = l + gammaln(1 + sum(x)) - sum(gammaln(1+x));
%    end;
%    l = l + sum(gammaln(p)) - gammaln(sum(p));
%  end;%
