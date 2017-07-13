function XC = center(X)
%CENTER XC = center(X)
%	Centers each column of X.


[T n] = size(X);
XC = X - ones(T,1)*(sum(X)/T); % Much faster than MEAN with a FOR loop
