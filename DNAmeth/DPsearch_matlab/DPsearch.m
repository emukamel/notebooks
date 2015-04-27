function [r,marginalOrder] = DPsearch(X,alpha,G0alpha,beamSize,Y,heuristic)
  % X         data set
  % X(n,w) = # of times w appears in trial n
  % alpha     concentration parameter for DP
  % G0alpha   concentration parameter for G0 prior
  % beamSize  size of the beam
  % Y         (optional) True clusters, defaults to ones(1,N)
  % heuristic 'n' no heuristic, 'a' admissible heuristic, 'i' inadmissible heuristic
    
[N,W] = size(X);

% initialize marginals

marginals = zeros(N,1);
for n=1:N,
  marginals(n) = log_marginal_posterior(X(n,:),G0alpha,[]);
end;

[X,Y,marginals,marginalOrder] = order_by_marginal(X,Y,marginals);

for n=1:(N+1),
  marginals(n) = sum(marginals(n:N));
end;

tic;
% initialize
s0.c        = 1;
s0.m        = [1];
s0.K        = 1;
s0.phi(1,:) = X(1,:);
s0.cnt(1)   = 1;
[s0.g,s0.l] = compute_scor(s0,X,alpha,G0alpha);
s0.h        = compute_heur(s0,X,alpha,G0alpha,marginals,heuristic);
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
  
	if N0 > bigK && mod(N0, 500) == 0,
    tm = toc;
    %[fs pr re] = compute_f(Y(1:N0), s.c);
    fs = 0; pr = 0; re = 0;
    disp(sprintf('%g \tn=%d k=%d |h|=%d f=%g', tm, N0, s.K, h.count,fs));
    bigK = N0;
  end;

  % check for completion
  if N0 == N,
    r = s;
    r.time = toc;
    r.numDequeued = numDequeued;
    r.numEnqueued = numEnqueued;
    %r.beam = h;
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
    [s1.g,s1.l] = compute_scor(s1,X,alpha,G0alpha);
    s1.h        = compute_heur(s1,X,alpha,G0alpha,marginals,heuristic);
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
  [s1.g,s1.l] = compute_scor(s1,X,alpha,G0alpha);
  s1.h        = compute_heur(s1,X,alpha,G0alpha,marginals,heuristic);
  s1.score    = s1.g + s1.h;
  s1.par      = [];
    
  h = heapinsert(h,s1); numEnqueued = numEnqueued+1;
  
  if h.count > beamSize,
    h.count = beamSize;
    h.tree  = h.tree(1:beamSize);
  end;
end;


function [l,dl] = compute_scor(s,X,alpha,G0alpha)
  N0 = length(s.c);
  if N0 == 1,
    dl = data_likelihood(X,G0alpha,s.c);
  else
    dl = s.par.l + data_likelihood_update(X,G0alpha,s.c,s.phi,s.cnt);
  end;    
  l  = - dl - log_DP_prior_count_complete2(s.c, alpha, size(X,1), s.m);

function l = compute_heur(s,X,alpha,G0alpha,marginals,heuristic)
  if heuristic == 'i',
    l = compute_heur_inad(s,X,alpha,G0alpha,marginals);
  elseif heuristic == 'a'
    l = compute_heur_admi(s,X,alpha,G0alpha,marginals);
  elseif heuristic == 'n',
    l = 0;
  end;
 

% 33.9519

function l = compute_heur_admi(s,X,alpha,G0alpha,marginals)
  [N,W] = size(X);
  N0    = length(s.c);
  K0    = s.K;

  badc = [9];
  
  lp = 0;
  if N0 < N,
    for n=(N0+1):N,
      x    = X(n,:);
      lik  = -Inf * ones(1, K0+1);
      
      for k=1:K0+1,
        post = X(find(s.c==k),:);
        lik(k) = log_marginal_posterior(x, G0alpha, post);
      end;
      [maxLik, k_max] = max(lik);
      for k=1:(K0+1), % [k_max (K0+1)],
        kept = [];
        post = X(find(s.c==k),:);
        for m=(N0+1):N,
          if m ~= n,
            if length(s.c) == length(badc) && all(s.c == badc) && n == 7,
              m
              [lik(k) lik2]
            end;
            lik2 = log_marginal_posterior(x, G0alpha, [post ; X(m,:)]);
            if lik2 > lik(k),
              lik(k) = lik2;
              kept = [kept m];
              post = [post ; X(m,:)];
            end;
          end;
        end;
      end;
      if length(s.c) == length(badc) && all(s.c == badc),
        n
        lik
        max(lik)
        kept
      end;
      lp = lp + max(lik);
    end;
  end;
  
  kv = log_DP_prior_count_complete2(s.c, alpha, N);

  if length(s.c) == length(badc) && all(s.c == badc),
    s.l
    lp
    kv
    s.g
%   fail;
  end;

  
  l = - s.l - lp - kv - s.g;
  
function l = compute_heur_inad(s,X,alpha,G0alpha,marginals)
  N0 = length(s.c);
  l  = - marginals(N0+1);

  %   - sl - m - kv - sg
  % = - sl - m - kv - (- dl - dp)
  % = - sl - m - kv + dl + dp
  % = - sl - m + dp