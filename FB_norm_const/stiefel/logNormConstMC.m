function [logC,se_logC] = logNormConstMC(A,B,n)

% LOGNORMCONSTMC  Log of Monte Carlo estimate of the normalising constant  
% of the Fisher-Bingham distribution on the Stiefel manifold, based on MC 
% samples drawn from the uniform distribution on the Stiefel manifold.
%
% A, B are parameters of the distribution (see the help text of
% logNormConstSP for details); n is the number of MC samples to use
% [default = 1e5].
%
% [logC,se_logC] = logNormConstMC(A,B,n) returns an estimate se_logC of 
% the standard error of logC.
%
% Example:
%    d=3; p=3; 
%    A = diag([10 5 1]); 
%    B = magic(d*p)/30; 
%    logC = logNormConstMC(A,B);
%
% See also: logNormConstSP

[d,p] = size(A);

if d<p, error('p cannot be larger than d'); end
if ~all(size(B) == d*p), error('B should be (d*p)-by-(d*p)'); end

if nargin<3, n=1e5; end

sample = zeros(n,1);

for k = 1:n
   Y = randOrthMat(d);
   X = Y(:,1:p);
   sample(k) = exp(trace(A(:)'*X(:))+ X(:)'*B*X(:));
end

% Multiply by the volume of the manifold
sample = sample * volumeOfStiefelManifold(d,p);

logC = log(mean(sample));

% Delta-method-based standard error
se_logC = sqrt(var(sample)/n)/mean(sample);

