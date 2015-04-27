from numpy import ones,zeros,argsort,asarray,floor,arange
from scipy.special import gammaln
import numpy as np

def DPsearch(X,alpha,G0alpha,beamSize,y=None,heuristic='i'):
	"""
	Arguments:
		X         matrix holding data (NxW), N trials, dimension = W
		alpha     concentration parameter for DP
		G0		  dict of parameters for beta prior distribution
				  G0.a = alpha parameters (vector, size W)
			      G0.b = beta parameters (vector, size W)
		beamSize  size of the beam (how many solutions to keep in queue)

	Optional arguments:
		Y         True clusters, defaults to ones(1,N)
		heuristic 'n' no heuristic, 'a' admissible, 'i' inadmissible
    """

    # Check inputs
    N,W = X.shape
    if y is none: y = ones(N)

    # Initialize marginal distributions:
    marginals = []
    for n in range(N):
    	marginals.append(beta_binomial_log_likelihood(X[n,:],G0))

    # Re-order data by marginal probabilities
    X,y,marginals = order_by_marginal(X,y,asarray(marginals))

    # marginals monotone increasing (why? and isn't this decreasing?)
    for n in range(N+1):
    	maginals[n] = np.sum(marginals[n:])

    # Initialize with one cluster, containing the first observation
    # 	n0 = 1 (number of datapoints from dataset we are considering)
    # 	each cluster defines a probability vector p[c,:]
    #   Each datapoint is modeled as X[n,i] ~ Bernoulli(p[j,i])
    #   m = vector containing number of clusters of each size.
    #   m[0],m[1],m[2] are number of clusters containing 1,2,3 datapoints
    #   N = # observations = sum(m[i]*(i+1))
    #   K = # clusters = sum(m[i]); 
    s0.c   = [1];      # cluster identities for first n0 datapoints 
    s0.m   = [1];      # vector of params/clusters drawn from DP (see above)
    s0.K   = 1;        # number of clusters?      
    s0.phi = X[0,:];   # PHI is sum(X[c,:],axis=0) in each cluster c
    s0.cnt = [1];      # Number of observations in each cluster
    
    # compute liklihood and heuristic for the initialized cluster
    # s0.l is log(H(x_1)), s0.g is -log[H(x_0|c_0)*max_{c0,c1}(p(c_0,c1))]
    # in inadmissible chace, -log \sum_{n=N0+1}^N H(x_n)
    [s0.g,s0.l] = compute_scor_betabernoulli(s0,X,alpha,G0)
    s0.h        = compute_heur_betabernoulli(s0,marginals,heuristic); 
    s0.score    = s0.g + s0.h;

def compute_scor_betabernoulli(s,X,alpha,G0):
	if len(s.c) == 1:
		# one cluster, directly compute data likelihood
		dl = beta_binomial_log_likelihood(X,G0);
	else:
		# data likelihood is parent likelihood plus and update term
		dl = s.par.l + data_likelihood_update_betabernoulli(X,G0,s.c,s.phi,s.cnt);

	l  = - dl - log_DP_prior_count_complete2(s.c, alpha, size(X,1), s.m); 
	return (l,dl)


def log_DP_prior_count_complete2(c0,alpha,N,m=None):
	if "dd1" not in log_DP_prior_count_complete2.__dict__:
		dd1  = arange(N-1)/(1+ arange(N-1))

		foo.counter = 0

	N0 = length(c0)
	if m is None: m = counts(c0,N)

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

def counts(c,N):
	"""
	this assumes that there are no "empty" clusters ...
	safe for search, unsafe for sampling (unless you GC every iteration)
	"""
	k = np.max(c) # number of clusters
	m = zeros(N)  
	t = zeros(k)
	for i in range(len(c)):
		t[c[i]] = t[c[i]] + 1
	for i in range(k):
		m[t[i]] = m[t[i]] + 1
	return m


def data_likelihood_update_betabernoulli(X,G0,c,phi,cnt):
  %[,W] = size(X);
  M     = length(c);
  
  if M <= 1,
    l = data_likelihood_betabernoulli(X,G0alphabeta);
  else
    k = c(M);
    
    if cnt(k) == 1, % new cluster
      l = data_likelihood_betabernoulli(X(M,:),G0alphabeta,1);
    else            % existing cluster
      xm = X(M,:);
      sx = phi(k,:);
      
     l = sum(gammaln(sx + G0alphabeta(1,:)) + gammaln(cnt(k)-sx+G0alphabeta(2,:))-...
         gammaln(sx-xm + G0alphabeta(1,:))- gammaln(cnt(k)-1 -(sx-xm)+G0alphabeta(2,:))+...
         gammaln(sum(G0alphabeta,1)+cnt(k)-1)-gammaln(sum(G0alphabeta,1)+cnt(k)));

def order_by_marginal(X,y,m):
	"""
	X 	data matrix, shape (NxW)
	Y 	clusters, shape (Nx1)
	M 	marginal probabilities for each observation, shape (Nx1)
	"""
	ind = argsort(m) # returns indices to sort ascending
	newX = X[ind,:]
	newY = y[ind]
	newM = m[ind]
	return (newX,newY,newM)

def beta_binomial_log_likelihood(X,G):
	"""
		X 	data (each row is an observation)
		G	parameters for beta distribution
	"""
	n_s,w_s = X.shape
	x_sum = np.sum(X, axis=0)
	return np.sum(gammaln(n_s + 1) + \
		    gammaln(x_sum + G.a) + \
	 		gammaln(n_s - x_sum + G.b) + \
	 		gammaln(G.a + G.b) - \
	 		gammaln(x_sum + 1) - \
	 		gammaln(n_s - x_sum + 1) - \
	 		gammaln(n_s + G.a + G.b) - \
	 		gammaln(G.b) - \
	 		gammaln(G.a))

def heap_insert(h,val):
	"""
		h 		heap structure (python dict), if h is None heap is initialized
		val		value to insert into the heap
	"""
	# initialize heap if necessary
	if h is None:
		h['count'] = 0
		h['tree'] = []

	# add new element
	h['tree'].append(val)
	h['count'] += 1

	# swap elements around to maintain heap structure
	cur,parent = h.count,floor(cur/2)
	found = False

	while not found:
		if h['tree'][parent].score < hout['tree'][cur].score:
			found = True
		else: # swap parent and current node
			parent_copy = h['tree'][parent]
			h['tree'][parent] = h['tree'][cur]
			h['tree'][cur] = parent_copy
			cur = parent
			parent = floor(cur/2)
		# check if we've reached the top of the heap
		if cur == parent:
			found = True

	# heap is now sorted correctly
	return h

def log_marginal_posterior(X,alpha,xs):
	N,W = X.shape
	l = data_likelihood(x, alpha, 1);

	return l

def order_by_marginal():
 	pass 
 	[newX,newY,newM,J] = order_by_marginal(X,Y,M)
  
  [newM,I] = sort(M,1,'ascend');
   newX = X(I,:);
   newY = Y(I);
   for i=1:length(I),
     J(I(i)) = i;
   end;

def compute_it(m,N0,N,dd1,logs,alpha):
""" WHAT IS THIS??? """

function l = compute_it(m,N0,N,dd1,logs,alpha)
  finishUp = 0;
  lastN    = 0;
  lM       = length(m);
  if m(lM) ~= 0
    m(lM+1) = 0;
    lM = lM + 1;
  end
    
  for n=(N0+1):N
      
    scores = dd1(1:lM-1) .* m(1:lM-1) ./ (m(2:lM)+1); % zhu: this is l/(l+1)*m_l/(m_{l+1}+1), the factor when n belongs to one of the existing cluster. 
    [val,idx] = max(scores); %zhu: idx is the number of datapoints an existing cluster have which are most possible for n to fall in; m(idx) is the number of clusters with idx datapoints.
    if val < alpha / (m(1)+1) % zhu: alpha/(m(1)+1), the factor when n belongs to a new cluster of its own.
      idx = 0;  %zhu: in this case, idx=0, indicates that n forms a new cluster.
    end

    m(idx+1) = m(idx+1) + 1; %zhu: If x_n falls in a existing cluster with idx points, m_{l+1} need to increase by one, and m_l need to decrease by one. 
    if idx > 0 %zhu: m_l needs to decrease by one.
      m(idx) = m(idx) - 1;
    end
    if lM == idx+1
        m(idx+2) = 0; 
        lM = idx+2; 
    end

    j = idx+1; %zhu: the l values.
    v = j/(j+1) * m(j)/(m(j+1)+1); %zhu: the factor l/(l+1)*m_l/(m_{l+1}+1)
    if j > 1 && v > alpha / (m(1) + 1) && m(j-1) == 0 && all(scores <= v) % zhu: whether the largest cluster is sufficiently large so we can do the simple grow.
          if all(m((j+1):end) == 0),
            if all((dd1(j-1:lM-1) .* m(j-1:lM-1) ./ (m(j:lM)+1)) <= v),
              finishUp = j;
              lastN = n;
              break;
            end;
          end;
    end
    
  end
  
  if finishUp > 0
%    fprintf(1, '%d ', lastN-N0);
    m(finishUp) = 0;
    m(finishUp + N - lastN) = 1;
  end
      
  m2 = find(m>0);
%  if max(m2) > length(dd1) || max(m2) > length(logs)
%    size(dd1)
%    size(logs)
%    m2
%  end
  l = sum(m(m2)) * log(alpha) - sum(m(m2) .* logs(m2)) - sum(gammaln(m(m2)+1)); 
  %zhu: this is P(m|alpha,N) without the constant part.
