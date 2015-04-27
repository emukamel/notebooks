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
