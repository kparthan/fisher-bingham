function out = subm(C,d,i,j)

% returns the (i,j)th sub-block of C, where d is the (linear) size of 
% each block, assumed square and equal

out = C(1+(i-1)*d:i*d, 1+(j-1)*d:j*d);