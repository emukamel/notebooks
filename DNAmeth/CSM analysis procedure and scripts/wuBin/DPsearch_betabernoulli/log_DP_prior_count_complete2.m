function l = log_DP_prior_count_complete2(c0,alpha,N,m)
  N0   = length(c0);
  if nargin < 4,
    m    = counts(c0,N);
  end;
  % maxD = find(m>0,1,'last');

  persistent dd1;
  persistent logs;
  persistent hashtable;
  if isempty(dd1) || length(logs) ~= N
    dd1  = (1:(N-1)) ./ (2:N);
    logs = log(1:N);
    hashtable = [];
  end

  h = hash(m,N);
  if h > length(hashtable) || hashtable(h) == 0
    hashtable(h) = compute_it(m,N0,N,dd1,logs,alpha);
  end
  l = hashtable(h);
  
  
 % zhu2012: line 25-39 was commented out because they are not called. 
% function m = counts_safe(c,N)
%   m = zeros(1,N);
%   k = unique(c);
%   for i=1:length(k),
%     j = sum(c==k(i));
%     m(j) = m(j) + 1;
%   end;
% 
% 
%   
% function v = log_rising(a,n)
%   v = 0;
%   for i=1:n,
%     v = v + log(a+i-1);
%   end;


  
