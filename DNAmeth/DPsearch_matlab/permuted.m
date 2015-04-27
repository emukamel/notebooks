function z = permuted(F)
  z = F;
  for f=1:size(F,1),
    r = ceil(f + rand * (size(F,1)-f));
    t = z(f,:);
    z(f,:) = z(r,:);
    z(r,:) = t;
  end;
