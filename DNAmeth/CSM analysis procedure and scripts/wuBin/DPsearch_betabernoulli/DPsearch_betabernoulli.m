function [r,marginalOrder] = DPsearch_betabernoulli(X,alpha,G0alphabeta,beamSize,Y,heuristic)
  % see DPgibbs
  % zhu note:
  % X:           is the data matrix (N by w), each row is one read, each X(n,i) is a binary
  %              value indicating whether site i in read n is methylated or not. 
  % alpha:       concentration parameter of the DP. 
  % G0alphabeta: a matrix of size 2 by w, the beta prior parameters for each
  %              p_j, j=1,...,w. The first row contains the alpha
  %              parameters and the second row contains the beta
  %              parameters.
  % beamSize:    a positive integeter. (How many solutions to retain in the
  %              queue).
  % Y:           true clusers if known. Otherwise, use ones(1,N). N is num. of data.
  % heuristic:   'n' if no heuristic, 'a' if using the admissible heuristic, 
  %             'i' is using the inadmissible heuristic. We will always use 'i'.              
      
N = size(X, 1);

% initialize marginals
marginals = zeros(N,1);
for n = 1 : N
  marginals(n) = log_marginal_posterior_betabernoulli(X(n,:),G0alphabeta,[]); %zhu: logH(x_n), n=1,...,N.
end

[X,Y,marginals,marginalOrder] = order_by_marginal(X,Y,marginals);%zhu: order logH(x_n) in ascending order.
%Zhu: note all results are corresponding to the ordered data, not the
%original data.
for n=1:(N+1)
  marginals(n) = sum(marginals(n:N)); % zhu: marginals is monotone increasing now.
end

tic;
% initialize
s0.c        = 1;
s0.m        = 1;
s0.K        = 1;
s0.phi(1,:) = X(1,:); %zhu: the ith row of phi is the sum of X in cluster i.
s0.cnt(1)   = 1; %zhu: the ith component is the size of cluster i.
[s0.g,s0.l] = compute_scor_betabernoulli(s0,X,alpha,G0alphabeta); % zhu: s0.l is log(H(x_1)), s0.g is -log[H(x_0|c_0)*max_{c0,1}(p(c_0,c1))]
s0.h        = compute_heur_betabernoulli(s0,marginals,heuristic);%zhu: in inadmissible chace, -log \sum_{n=N0+1}^N H(x_n)
s0.score    = s0.g + s0.h;

h = heapinit;
h = heapinsert(h,s0);

bigK = 0;
% iterate
numDequeued = 0; numEnqueued = 0;
while (1),
  [s,h] = heapremove(h);
  numDequeued = numDequeued+1;
  N0 = length(s.c);
  
  if N0 > bigK && mod(N0, 500) == 0
    tm = toc;
    %[fs pr re] = compute_f(Y(1:N0), s.c);
    fs = 0; pr = 0; re = 0;
    %disp(sprintf('%g \tn=%d k=%d |h|=%d f=%g', tm, N0, s.K, h.count,fs));
    bigK = N0;
  end

  % check for completion
  if N0 == N,
    r = s;
    r.time = toc;
    r.numDequeued = numDequeued;
    r.numEnqueued = numEnqueued;
%    r.beam = h;
    [fs pr re] = compute_f(Y, s.c);
    r.prf = [fs pr re];
    break;
  end;
  
  % expand by existing cluster
  for k=1:s.K,
    sze = sum(s.c == k);
    m   = s.m;
    m(sze) = m(sze) - 1;
    if sze == length(m), m = [m 0]; end;
    m(sze+1) = m(sze+1) + 1;

    phi = s.phi; phi(k,:) = phi(k,:) + X(N0+1,:);
    cnt = s.cnt; cnt(k) = cnt(k)+1;
    
    s1.c        = [s.c k];
    s1.m        = m;
    s1.K        = s.K;
    s1.phi      = phi;
    s1.cnt      = cnt;
    s1.par      = s;
    [s1.g,s1.l] = compute_scor_betabernoulli(s1,X,alpha,G0alphabeta);
    s1.h        = compute_heur_betabernoulli(s1,marginals,heuristic);
    s1.score    = s1.g + s1.h;
    s1.par      = [];
    
    h = heapinsert(h,s1); numEnqueued = numEnqueued+1;
  end;
  
  % expand by new cluster
  phi = s.phi; phi(s.K+1,:) = X(N0+1,:);
  cnt = [s.cnt 1];
  
  s1.c        = [s.c (s.K+1)];
  s1.m        = s.m; s1.m(1) = s1.m(1) + 1;
  s1.K        = s.K+1;
  s1.phi      = phi;
  s1.cnt      = cnt;
  s1.par      = s;
  [s1.g,s1.l] = compute_scor_betabernoulli(s1,X,alpha,G0alphabeta);
  s1.h        = compute_heur_betabernoulli(s1,marginals,heuristic);
  s1.score    = s1.g + s1.h;
  s1.par      = [];
    
  h = heapinsert(h,s1); numEnqueued = numEnqueued+1;
  
  if h.count > beamSize,
    h.count = beamSize;
    h.tree  = h.tree(1:beamSize);
  end;
end;



 

 
  %   - sl - m - kv - sg
  % = - sl - m - kv - (- dl - dp)
  % = - sl - m - kv + dl + dp
  % = - sl - m + dp