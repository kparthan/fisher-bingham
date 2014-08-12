function out = subv(v,dcum,i)

% returns the (i)th sub-vector of v, where 
%     v = [va{1}, ... , va{p}];
% and
%     va{i} is da(i)-by-1
% where
%    d = [da(1); ... ; da(p)];
% and 
%    dcum = [0 cumsum(da)];

out = v(1+dcum(i):dcum(i+1));