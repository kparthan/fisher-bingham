function X = randVonMisesFisher(d,lambda,n)

% RANDVONMISESFISHER  Returns a vector drawn from the von Mises-Fisher
% distribution on the unit sphere S^{d-1} in R^{d}, with modal direction 
% (0,...,1)' and concentration parameter lambda (>= 0).  Parameter n is 
% the number of draws [default = 1].
%
% This is an implementation of the rejection-sampling algorithm in
% Wood (1994), "Simulation of the von Mises-Fisher Distribution", Commun.
% Statist. Simula. 23(1):157-164.

% Written by Simon Preston (http://www.maths.nott.ac.uk/~spp), 2010

if nargin<3, n = 1; end
if lambda<0, error('lambda must be non-negative'); end

if n>1
   X = zeros(d,n);
   for n_ = 1:n
      X(:,n_) = randVonMisesFisher(d,lambda,1);
   end
   return
end

b = (-2*lambda + sqrt(4*lambda^2+(d-1)^2))/(d-1);
x0 = (1-b)/(1+b);
c = lambda*x0 + (d-1)*log(1-x0^2);

rejCond = -1;
while rejCond < 0
   Z = betarnd((d-1)/2,(d-1)/2);
   U = rand;
   W = (1-(1+b)*Z)/(1-(1-b)*Z);
   rejCond = lambda*W + (d-1)*log(1-x0*W) - c - log(U);
end

V = randn(d-1,1);
V = V./norm(V);

X = [sqrt(1-W^2)*V; W];
