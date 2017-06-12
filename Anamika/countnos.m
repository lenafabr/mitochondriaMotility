function [outarray] = countnos(L,inarray)
%Function to count # of occurences in a given interval [l,l+1] for all
%intervals
% To get a distribution of mitochondria #s in the length [0,L]
%outarray : gives #s for all l in [0,L]
%Similar to histogram function, but with boundaries [0,L] irrespective of
%values, and fixed bin size of 1.
outarray = zeros(1,L); 
y = round(inarray);
for i = 1:1:size(y)
    ind = y(i);
    outarray(ind) = outarray(ind) + 1;
end
end

