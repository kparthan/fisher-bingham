function out = subm2(C,dcum,i,j)

% returns the (i,j)th sub-block of C, i.e., let Ca{i,j} be da(i)-by-da(j) 
%     C = [Ca{1,1}, ... , C{a1,p}; ...
%                   ...          
%          Ca{p,1}, ... , C{ap,p}];
% and
%    d = [da(1); ... ; da(p)];
% and 
%    dcum = [0 cumsum(d)];

out = C(1+dcum(i):dcum(i+1), 1+dcum(j):dcum(j+1));