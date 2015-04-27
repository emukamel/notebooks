function v = hash(m,N)
  persistent logs; %zhu2012: logs is persistent for this function, not for other functions.
  if isempty(logs) || length(logs) ~= N
    logs = log(7 + 3 * (1:N));
  end
  
  v = mod(floor(sum(m .* logs(1:length(m)))), 97) + 1;