function out = logNormConstSP(A,B,approxType)

% LOGNORMCONSTSP  Computes the log of the saddlepoint approximation to the 
% normalising constant of the Fisher-Bingham distribution on the Stiefel 
% manifold V_{d,p}, defined as the set {X : X'*X = eye(p)} where X is 
% d-by-p.
%
% Inputs:
%
%    - A, a d-by-p matrix of parameters.
%
%    - B, a (d*p)-by-(d*p) matrix of parameters.
%
%    - approxType: used for specifying the output to be
%           (1) first-order SP, (2) second-order SP "1+T" variant,
%           (3) second-order SP "exp(T)" variant, or (4) a vector
%           containing all three of the above [default = 3]
%
% Example:
%    d=3; p=3; 
%    A = diag([10 5 1]); 
%    B = magic(d*p)/30; 
%    logC = logNormConstSP(A,B);
% 

% Reference:
%    Kume, A., Preston, S.P., Wood, A.T.A., "Saddlepoint approximations for
%    the normalising constant of Fisher-Bingham distributions on products
%    of spheres and Stiefel manifolds", draft manuscript, 2011

if nargin < 3 || isempty(approxType), approxType = 3; end

[d,p] = size(A);

if d<p, error('p cannot be larger than d'); end
if ~all(size(B) == d*p),
   error('B should be (d*p)-by-(d*p)');
end

% make sure that B is symmetric and "-ve def" (so that V is "+ve def")
B = 0.5*(B+B');
evals = eig(B);
maxRealPart = max(real(evals));
if maxRealPart >= 0
   B = B - 2*(maxRealPart+1)*eye(d*p);
   logCorrectionFactor = 2*(maxRealPart+1)*p;
else
   logCorrectionFactor = 0;
end

V = inv(-2*B);
mu = (-2*B)\vec(A); %V*vec(A);

theta_v_initialGuess = upperTri(zeros(p));

optimOptions = optimset('Display','off','TolX',1e-30,'TolFun',1e-30,'LargeScale','off');

thetaHat_v = fminunc(@(theta_v) objFunc(theta_v,p,d,mu,V), ...
   theta_v_initialGuess, optimOptions);

thetaHat_m = upperTri(thetaHat_v);

% % Check output of optimisation
% display(' '); display(' ');
% display('Check K1 - eye(p) is close to zero:');
% display(upperTri(calc_K1(thetaHat_m,p,d,mu,V)) - eye(p));

% lambda_m = V\mu;
CHat = calc_C(thetaHat_m,p,d,V);
CHatInv = inv(CHat);

Khat = calc_K(thetaHat_m,p,d,mu,V);
K2hat = calc_K2(thetaHat_m,p,d,mu,V);

%K1hat = calc_K1(thetaHat,p,d,mu,V); % not needed - only for checking purposes

m = 1/2*p*(p+1);

log_f2Hat = -m/2*log(2*pi) - 1/2*log(det(K2hat)) + Khat - sum(diag(thetaHat_m));

log_C2_firstOrder = logCorrectionFactor + log_f2Hat + (d*p/2)*log(2*pi) + ...
   1/2*log(det(V)) + 1/2 * mu'/V*mu + p*log(2);

% log_C2_firstOrderVersion2 = logCorrectionFactor + p*log(2) + ...
%    (d*p/2 - p*(p+1)/4)*log(2*pi) - 1/2*log(det(K2hat)) - 1/2*log(det(CHat)) + ...
%    1/2*(vec(A)')/CHat*vec(A) - sum(diag(thetaHat_m));

if approxType== 1
   out = log_C2_firstOrder;
else
   T = calc_T(CHatInv,p,d,mu,V);
   switch approxType
      case 2, out = log_C2_firstOrder + log(1+T);
      case 3, out = log_C2_firstOrder + T;
      case 4, out = [log_C2_firstOrder,...
            log_C2_firstOrder + log(1+T), ...
            log_C2_firstOrder + T];
   end
end


% ---------------------------------------------------------------------

function [out1, out2] = objFunc(theta_v,p,d,mu,V)

theta = upperTri(theta_v);

C = calc_C(theta,p,d,V);
[evecs,D] = eig(C);
if min(diag(D)) < 1e-14
   out1 = 1e99;
else
   out1 = calc_K(theta,p,d,mu,V) - sum(diag(theta));
end

