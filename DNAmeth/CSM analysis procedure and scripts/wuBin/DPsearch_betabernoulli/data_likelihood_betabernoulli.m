function l = data_likelihood_betabernoulli(X,G0alphabeta,c)
  %[~, W] = size(X);
  if length(c) > 1,
    %fail; zhunote: for this case, the author did not code. 
    % the following code is added by zhu.
    x_sum = sum(X, 1);
    n_s = size(X,1);
%     l = gammaln(W*alpha) - W*gammaln(alpha) + ...
%         sum(gammaln(alpha+x_sum)) - gammaln(sum(x_sum+alpha)) + ...
%         sum(gammaln(sum(X,2)+1)) - sum(sum(gammaln(X+1)));
      l = sum(gammaln(x_sum+G0alphabeta(1,:)) + gammaln(n_s - x_sum + G0alphabeta(2,:)) - gammaln(G0alphabeta(1,:)) - gammaln(G0alphabeta(2,:)) ...
          + gammaln(sum(G0alphabeta,1)) - gammaln(sum(G0alphabeta,1)+n_s));      
  elseif length(c) == 1
  
%     l = gammaln(W*alpha) - W*gammaln(alpha) + ...
%         sum(gammaln(alpha + X(1,:))) - gammaln(sum(alpha + X(1,:))) + ...
%         gammaln(sum(X(1,:)) + 1) - sum(gammaln(X(1,:)+1));
      l = sum(gammaln(X(1,:)+G0alphabeta(1,:)) + gammaln(1-X(1,:)+G0alphabeta(2,:))- gammaln(G0alphabeta(1,:)) - gammaln(G0alphabeta(2,:)) +...
          gammaln(sum(G0alphabeta,1)) - gammaln(G0alphabeta(1,:)+G0alphabeta(2,:)+1));
  end
    
  
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
