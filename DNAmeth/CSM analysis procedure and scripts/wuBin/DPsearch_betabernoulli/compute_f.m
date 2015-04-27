function [f,p,r] = compute_f(T,H)

% Zhu:
% f: the f-score. The geometric mean of precision and recall. 
% p: precision. How many of element of the cluster in H are truely a
%               cluster.
% r: recall. did all the element cluster in T make it.
% reference: Zhao and Karypis, 2002. Criterion functions for document
% clustering. Technical report, experiments and analysis University of
% Minnesota, Department of Computer Science/Army HPC research center. 
  if length(T) ~= length(H),
    size(T)
    size(H)
  end;
  
  N = length(T);
  numT = 0;
  numH = 0;
  numI = 0;
  for n=1:N,
    Tn = (T(n+1:end))==T(n); %zhu2012:  not sure what this is doing.
    Hn = (H(n+1:end))==H(n);
    numT = numT + sum(Tn);
    numH = numH + sum(Hn);
    numI = numI + sum(Tn .* Hn);
%    for m=(n+1):N,
%      if T(n) == T(m),
%        numT = numT+1;
%        if H(n) == H(m),
%          numI = numI+1;
%        end;
%      end;
%      if H(n) == H(m),
%        numH = numH+1;
%      end;
%    end;
  end;
  p = 1;
  r = 1;
  f = 1;
  if numH > 0,
    p = numI / numH;
  end;
  if numT > 0,
    r = numI / numT;
  end;
  if (p+r) == 0,
    f = 0;
  else
    f = 2 * p * r / (p + r);
  end;
  