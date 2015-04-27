function l = compute_heur_betabernoulli(s,marginals,heuristic)
  if heuristic == 'i',
    l = compute_heur_inad_betabernoulli(s,marginals);
  elseif heuristic == 'a'
   % l = compute_heur_admi(s,X,alpha,G0alpha,marginals);
  elseif heuristic == 'n',
   % l = 0;
  end;
