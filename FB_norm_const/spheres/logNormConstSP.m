function logC = logNormConstSP(d,a,B,approxType)

% LOGNORMCONSTSP  Computes the log of the saddlepoint approximation to the 
% normalising constant of the Fisher-Bingham distribution on the product 
% of spheres.
%
% Inputs:
%
%    - d = [d{1} ; ... ; d{p}];  where d{i} is a positive integer
%
%    - a = [a{1} ; ... ; a{p}];  where a{i} is d{i}-by-1
%
%    - B = [B{1,1}, ... , B{1,p}; ...
%                   ...
%           B{p,1}, ... , B{p,p}]  where B{i,j} is d{i}-by-d{j}
%
%     - approxType: used for specifying the output to be
%           (1) first-order SP, (2) second-order SP "1+T" variant,
%           (3) second-order SP "exp(T)" variant, or (4) a vector
%           containing all three of the above [default = 3]
%
% Example:
%    d = [2 2]';
%    a = [10 0 5 0]';
%    B = zeros(4); B(2,4) = 1;
%    logC = logNormConstSP(d,a,B);

% Notes:
%    The 'cell-array'-style notation is used above only to schematically 
%    illustrate the block structure of the inputs. Inputs a and d are 
%    column vectors and B is a matrix.

% Reference:
%    Kume, A., Preston, S.P., Wood, A.T.A., "Saddlepoint approximations for
%    the normalising constant of Fisher-Bingham distributions on products
%    of spheres and Stiefel manifolds", draft manuscript, 2011

if nargin < 4 || isempty(approxType), approxType = 3; end

p = length(d);  % number of spheres
dcum = [0 cumsum(d')];
D = sum(d);

% ensure that B is symmetric and "-ve def" (so that V is "+ve def"); if 
% not, use the distribution's degeneracy property to adjust accordingly
B = 0.5*(B+B');
evals = eig(B);
maxRealPart = max(real(evals));
if maxRealPart == 0, maxRealPart = 0.1; end;
if maxRealPart >= 0
   B = B - 2*maxRealPart*eye(D);
   logCorrectionFactor = 2*maxRealPart*p;
else
   logCorrectionFactor = 0;
end

V = inv(-2*B);
mu = (-2*B)\a;  % i.e. mu = V*a

theta0 = zeros(p,1);  % starting point for optimisation

optsDerivs = optimset('Display','notify','GradObj','on','Hessian','off','TolX',1e-20,...
   'TolFun',1e-20,'MaxFunEvals',1e5,'MaxIter',1e5);

optsNoDerivs = optimset('Display','notify','GradObj','off','LargeScale','off','TolX',...
   1e-17,'TolFun',1e-17,'MaxFunEvals',1e5,'MaxIter',1e5,'Diagnostics','off');


% try to make use of analytic Jacobian (doing so in fminunc in 
% older versions of Matlab can fail due to non-robust handling of the NaNs 
% which are returned by the objective function for non-feaasible theta); 
% otherwise fall back on non-Jacobian version
try
   thetaHat = fminunc(@(theta) objFuncDerivs(theta,p,dcum,mu,V), theta0, optsDerivs);
catch  %#ok<CTCH>
   thetaHat = fminunc(@(theta) objFuncNoDerivs(theta,p,dcum,mu,V), theta0, optsNoDerivs);
end

CHat = calc_C(thetaHat,p,dcum,V);

% % Check output of optimisation
% display(' '); display(' ');
% display('Check:  K1 - 1 =');
% display(num2str(calc_K1(CHat,p,dcum,mu,V) - ones(p,1)));

log_f = -p/2*log(2*pi) -1/2*log(det(calc_K2(CHat,p,dcum,mu,V))) + ...
   calc_K(CHat,p,dcum,mu,V) - thetaHat'*ones(p,1);

log_C1_firstOrder = logCorrectionFactor + log_f + (sum(d)/2)*log(2*pi) + ...
   p*log(2) + 1/2*log(det(V)) + 1/2*mu'/V*mu;

if approxType == 1
   logC = log_C1_firstOrder;
else
   T = 1/8*calc_rho4(CHat,p,dcum,mu,V) - ...
      1/24*(3*calc_rho13sq(CHat,p,dcum,mu,V) + ...
      2*calc_rho23sq(CHat,p,dcum,mu,V));
   switch approxType
      case 2, logC = log_C1_firstOrder + log(1+T);
      case 3, logC = log_C1_firstOrder + T;
      case 4, logC = [log_C1_firstOrder,...
            log_C1_firstOrder + log(1+T), ...
            log_C1_firstOrder + T];
   end
end


% ---------------------------------------------------------------------

function out = objFuncNoDerivs(theta,p,dcum,mu,V)

C = calc_C(theta,p,dcum,V);
[foo,evals_m] = eig(C);
if min(diag(evals_m)) < 1e-14
   out = 1e99;
else
   Cinv = inv(C);
   out = calc_K(Cinv,p,dcum,mu,V) - sum(theta);
end


% ---------------------------------------------------------------------

function [out1, out2, out3] = objFuncDerivs(theta,p,dcum,mu,V)

C = calc_C(theta,p,dcum,V);

[foo,evals_m] = eig(C);
if min(diag(evals_m)) < 1e-14
   out1 = NaN;
   out2 = NaN;
   out3 = NaN;
else
   out1 = calc_K(C,p,dcum,mu,V) - sum(theta);
   if nargout > 1, out2 = calc_K1(C,p,dcum,mu,V) - 1; end
   if nargout > 2, out3 = calc_K2(C,p,dcum,mu,V); end
end
