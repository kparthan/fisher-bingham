function [r,s] = get_r_s(k)

% Suppose that k is the index of the vectorized elements of an upper
% triangular square matrix, e.g. of the form,
%    1 2 4 7
%    0 3 5 8
%    0 0 6 9
%    0 0 0 10
% then [r,s] are the corresponding co-ordinates.  
% E.g get_r_s(9) returns [3,4].

s = floor(sqrt(2*k - 7/4) + 1/2);
r = k - 1/2*s*(s-1);