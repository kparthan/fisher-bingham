function T = calc_T(CHatInv,p,d,mu,V,matrixFisherFlag)

if nargin < 6 || isempty(matrixFisherFlag), matrixFisherFlag = 0; end

lambda = V\mu;

if matrixFisherFlag
   foo = diag(CHatInv);
   phiHat_d = foo(1:d:end);
   lambda = V\mu;
   gam = lambda(1:d+1:end);
   thetaHat = diag(1/2*(1-2*gam.^2./(sqrt(d^2+4*gam.^2)-d)));
   K2HatInv = inv(calc_K2(thetaHat,p,d,mu,V));
   K3Hat = calc_KHat3_matrixFisher(phiHat_d,gam,d);
   K4Hat = calc_KHat4a_matrixFisher(phiHat_d,gam,d);
else
   K2HatInv = inv(calc_KHat_n(2,CHatInv,p,d,[],lambda));
   K3Hat = calc_KHat_n(3,CHatInv,p,d,[],lambda);
   K4Hat = calc_KHat_n(4,CHatInv,p,d,[],lambda);
end

m = 1/2*p*(p+1);

%%%% rho4 %%%%

rho4 = 0;

for i = 1:m
   for j = 1:m
      for k = 1:m
         for l = 1:m
            indSrt = sort([i j k l]);  % sorted indices
            rho4 = rho4 + ...
               K4Hat(indSrt(1),indSrt(2),indSrt(3),indSrt(4))*...
               K2HatInv(i,j)*K2HatInv(k,l);
         end
      end
   end
end

%%%% rho13sq %%%%

rho13sq = 0;

for i = 1:m
   for j = 1:m
      for k = 1:m
         for r = 1:m
            for s = 1:m
               for t = 1:m
                  indSrtA = sort([i j k]);
                  indSrtB = sort([r s t]);
                  rho13sq = rho13sq + ...
                     K3Hat(indSrtA(1),indSrtA(2),indSrtA(3)) * ...
                     K3Hat(indSrtB(1),indSrtB(2),indSrtB(3)) * ...
                     K2HatInv(i,j)*K2HatInv(k,r)*K2HatInv(s,t);
               end
            end
         end
      end
   end
end

%%%% rho23sq %%%%

rho23sq = 0;

for i = 1:m
   for j = 1:m
      for k = 1:m
         for r = 1:m
            for s = 1:m
               for t = 1:m
                  indSrtA = sort([i j k]);
                  indSrtB = sort([r s t]);
                  rho23sq = rho23sq + ...
                     K3Hat(indSrtA(1),indSrtA(2),indSrtA(3)) * ...
                     K3Hat(indSrtB(1),indSrtB(2),indSrtB(3)) * ...
                     K2HatInv(i,r)*K2HatInv(j,s)*K2HatInv(k,t);
               end
            end
         end
      end
   end
end

% check:
% display([rho4 rho13sq rho23sq]);

T = rho4/8 - 1/24*(3*rho13sq + 2*rho23sq);


