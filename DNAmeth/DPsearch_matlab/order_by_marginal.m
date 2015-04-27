function [newX,newY,newM,J] = order_by_marginal(X,Y,M)
  
  [newM,I] = sort(M,1,'ascend');
   newX = X(I,:);
   newY = Y(I);
   for i=1:length(I),
     J(I(i)) = i;
   end;
   