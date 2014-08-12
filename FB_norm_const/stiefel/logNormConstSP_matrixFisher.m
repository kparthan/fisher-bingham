function out = logNormConstSP_matrixFisher(A,approxType)

% LOGNORMCONSTSP_MATRIXFISHER  Computes the log of the saddlepoint 
% approximation to the normalising constant of the Fisher matrix 
% distribution on the Stiefel manifold, V_{d,p}, defined as the set 
% {X : X'*X = eye(p)} where X is d-by-p.  
%
% Inputs:
%
%    - A, a d-by-p matrix of parameters.
%
%    - approxType: used for specifying the output to be 
%           (1) first-order SP, (2) second-order SP variant 1, 
%           (3) second-order SP variant 2, or (4) a vector containing 
%           all three of the above [default = 3]

% Notes:  
%     This distribution is a special (B=0) case of the Fisher-Bingham 
%     distribution on the Stiefel manifold.  In this case the saddlepoint 
%     equation is solvable in closed form, i.e. requires no numerical 
%     optimisation.   Calculations for the second-order approximations 
%     are also much simpler and faster than for the more general 
%     Fisher-Bingham distribution; see the reference below for details.

% Reference:
%    Kume, A, Preston, S.P., Wood, A.T.A, "Saddlepoint approximations for
%    the normalising constant of Fisher-Bingham distributions on products
%    of spheres and Stiefel manifolds", draft manuscript, 2011

if nargin < 2 || isempty(approxType), approxType = 3; end

[d,p] = size(A);
[Q,lambda_m] = svd(A);
lambda_d = diag(lambda_m);

% Here we take V to be the identity.  This corresponds to taking
% B = -1/2 * eye(d*p) (see SPP notes 6/7/11).  The property satisfied by C
% is C(A,B+kron(D,eye(p))) = C(A,B)*exp(trace(D)), i.e.,
% C(A,0) = exp(-trace(D))*C(A,kron(D,eye(p))),
% so in the following calculations, D = -1/2*eye(p), hence the "correction
% factor" of exp(-trace(D)) = exp(p/2).
logCorrectionFactor = p/2;
V = eye(p*d);

% note: A1 is the version of A diagonalised by applying an isometry
A1 = zeros(size(A));
A1(1:p,1:p) = diag(lambda_d);
mu = A1(:);

phiHat_d = (-p + sqrt(p^2+4*lambda_d.^2))./(2*lambda_d.^2);
phiHat_d(isnan(phiHat_d)) = 0; % accounts for case in which lambda_d = 0

thetaHat_d = 1/2*(1 - 2*lambda_d.^2./(sqrt(p^2+4*lambda_d.^2)-p));
thetaHat = diag(thetaHat_d);
ThetaTilde_m = 2*thetaHat;

CHat = kron((eye(p)-ThetaTilde_m),eye(d));

Khat = -.5*log(det(CHat)) - .5*log(det(V)) + ...
    .5*(mu'/(V*CHat*V)*mu) -.5*(mu'/V*mu);

K2hat = calc_KHat2_matrixFisher(phiHat_d,lambda_d,d);

CHatInv = inv(calc_C(thetaHat,p,d,V));

m = 1/2*p*(p+1);
log_f2Hat = -m/2*log(2*pi) - 1/2*log(det(K2hat)) + Khat - sum(diag(thetaHat));

log_C2_firstOrder = logCorrectionFactor + log_f2Hat + (d*p/2)*log(2*pi) + ...
   1/2*log(det(V)) + 1/2 * mu'/V*mu + p*log(2);

out = log_C2_firstOrder;

if approxType== 1
   out = log_C2_firstOrder;
else
   T = calc_T_matrixFisher(CHatInv,p,d,mu,V);
   switch approxType
      case 2, out = log_C2_firstOrder + log(1+T);
      case 3, out = log_C2_firstOrder + T;
      case 4, out = [log_C2_firstOrder,...
            log_C2_firstOrder + log(1+T), ...
            log_C2_firstOrder + T];
   end
end