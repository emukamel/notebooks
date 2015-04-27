function [newX,newY,newM,J] = order_by_marginal(X,Y,M)
  
  [newM,I] = sort(M,1,'ascend'); % zhu: newM = M(I)
   newX = X(I,:);
   newY = Y(I);
%    for i=1:length(I),
%     J(I(i)) = i;
%    end;
%    
   % zhu2012: replace line 6-8 by 11-12. 
   J=NaN(1,length(I));
   J(I)=1:length(I); %zhu: J indicate the rank of each component of M.
   