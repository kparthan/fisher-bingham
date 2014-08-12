function [logC, se_logC] = logNormConstMC(d,a,B,n)

% LOGNORMCONSTMC  Log of Monte Carlo estimate of the normalising constant 
% of the Fisher-Bingham distribution on the product of spheres, based on MC
% samples drawn from the uniform distribution on the product of spheres.
%
% d, a, B are parameters of the distribution (see the help text of
% logNormConstSP for details); n is the number of MC samples to use
% [default = 1e5].
%
% [logC,se_logC] = logNormConstMC(A,B,n) returns an estimate se_logC of 
% the standard error of logC.
%
% Example:
%    d = [2 2]';
%    a = [10 0 5 0]';
%    B = zeros(4); B(2,4) = 1;
%    logC = logNormConstMC(d,a,B);
%
% See also: logNormConstSP

if nargin<4 || isempty(n), n = 1e5; end

p = length(d);
D = sum(d);
dcum = [0 cumsum(d')];

% random sample on product of spheres
X = zeros(D,n);
surfaceAreaOfSphere = zeros(p,1);
for k = 1:p
   X(1+dcum(k):dcum(k+1),:) = randUniformOnSphere(d(k),n);
   surfaceAreaOfSphere(k) = 2*pi^(d(k)/2)/gamma(d(k)/2);
end

% MC sample of the density function evaluated at the points defined in X
G = zeros(n,1);
for k = 1:n
   G(k) = exp(a'*X(:,k) + X(:,k)'*B*X(:,k));
end

% Multiply by the volume of the manifold
G = G * prod(surfaceAreaOfSphere);

logC = log(mean(G));

% Delta-method-based standard error
se_logC = sqrt(var(G)/n)/mean(G);
