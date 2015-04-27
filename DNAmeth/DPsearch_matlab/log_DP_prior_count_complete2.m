function l = log_DP_prior_count_complete2(c0,alpha,N,m)
  N0   = length(c0);
  if nargin < 4,
    m    = counts(c0,N);
  end;
  % maxD = find(m>0,1,'last');

  persistent dd1;
  persistent logs;
  persistent hashtable;
  if isempty(dd1) || length(logs) ~= N,
    dd1  = (1:(N-1)) ./ (2:N);
    logs = log(1:N);
    hashtable = [];
  end;

  h = hash(m,N);
  if h > length(hashtable) || hashtable(h) == 0,
    hashtable(h) = compute_it(m,N0,N,dd1,logs,alpha);
  end;
  l = hashtable(h);
  
  
function l = compute_it(m,N0,N,dd1,logs,alpha)
  finishUp = 0;
  lastN    = 0;
  lM       = length(m);
  if m(lM) ~= 0,
    m(lM+1) = 0;
    lM = lM + 1;
  end;
    
  for n=(N0+1):N,
    scores = dd1(1:lM-1) .* m(1:lM-1) ./ (m(2:lM)+1);
    [val,idx] = max(scores);
    if val < alpha / (m(1)+1),
      idx = 0;
    end;

    m(idx+1) = m(idx+1) + 1;
    if idx > 0,
      m(idx) = m(idx) - 1;
    end;
    if lM == idx+1, m(idx+2) = 0; lM = idx+2; end;

    j = idx+1;
    v = j/(j+1) * m(j)/(m(j+1)+1);
    if j > 1 && v > alpha / (m(1) + 1) && m(j-1) == 0 && all(scores <= v),
      if all(m((j+1):end) == 0),
        if all((dd1(j-1:lM-1) .* m(j-1:lM-1) ./ (m(j:lM)+1)) <= v),
          finishUp = j;
          lastN = n;
          break;
        end;
      end;
    end;
  end;
  
  if finishUp > 0,
%    fprintf(1, '%d ', lastN-N0);
    m(finishUp) = 0;
    m(finishUp + N - lastN) = 1;
  end;
      
  m2 = find(m>0);
  if max(m2) > length(dd1) || max(m2) > length(logs),
    size(dd1)
    size(logs)
    m2
  end;
  l = sum(m(m2)) * log(alpha) - sum(m(m2) .* logs(m2)) - sum(gammaln(m(m2)+1));
  
function m = counts_safe(c,N)
  m = zeros(1,N);
  k = unique(c);
  for i=1:length(k),
    j = sum(c==k(i));
    m(j) = m(j) + 1;
  end;

function m = counts(c,N)  % this assumes that there are no "empty"
                          % clusters ... safe for search, unsafe for
                          % sampling (unless you GC every iteration)
  k = max(c);
  m = zeros(1,N);
  t = zeros(1,k);
  for i=1:length(c), t(c(i)) = t(c(i)) + 1; end;
  for i=1:k, m(t(i)) = m(t(i)) + 1; end;
  
  
function v = log_rising(a,n)
  v = 0;
  for i=1:n,
    v = v + log(a+i-1);
  end;

function v = hash(m,N)
  persistent logs;
  if isempty(logs) || length(logs) ~= N,
    logs = log(7 + 3 * (1:N));
  end;
  
  v = mod(floor(sum(m .* logs(1:length(m)))), 97) + 1;
  
