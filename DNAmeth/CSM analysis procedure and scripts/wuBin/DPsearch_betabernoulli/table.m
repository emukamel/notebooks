function y=table(x)
% This function list frequency of vector x.
% The first row is the items in x, the second row contains the frequency
% counts.

names=unique(x);
y=NaN(2,length(names));

y(1,:)=names;
for i=1:length(names)
    y(2,i)=sum(x==names(i));
end
