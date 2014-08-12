function T = calc_T_matrixFisher(CHatInv,p,d,mu,V)

foo = diag(CHatInv);
phiHat_d = foo(1:d:end);

lambda = V\mu;
omega = lambda(1:d+1:end);
h = omega.^2.*phiHat_d/d;

sigma = zeros(7,1);


%%%% rho13sq

runsum = 0;
for r=1:p
   foo = 0;
   for s = 1:p
      if r ~= s
         foo = foo + (1+2*h(r)+h(s))/(1+h(r)+h(s));
      end
   end
   runsum = runsum + 1/(1+2*h(r))*(2*(1+3*h(r))/(1+2*h(r))+foo)^2;
end

rho13sq = runsum*2/d;


%%%% rho23sq

sigma(1) = sum((1+3.*h).^2./(1+2.*h).^3);

for r = 1:p
   for s = 1:p
      if s ~= r
         sigma(2) = sigma(2) + (1+2*h(r)+h(s))^2/((1+2*h(r))*(1+h(r)+h(s))^2);
      end
   end
end

for r = 1:p
   for s = r+1:p
      for t = s+1:p
         sigma(3) = sigma(3) + (1+h(r)+h(s)+h(t))^2/((1+h(r)+h(s))*(1+h(r)+h(t))*(1+h(s)+h(t)));
      end
   end
end

rho23sq = 2/d*(4*sigma(1)+3*sigma(2)+3*sigma(3));


%%%% rho23sq

sigma(4) = sum((1+4*h)./(1+2*h).^2);

for r = 1:p
   for s = r+1:p
      sigma(5) = sigma(5) + (1+2*h(r)+2*h(s))/(1+h(r)+h(s))^2;
   end
end

for r = 1:p
   for s = 1:p
      if s ~= r
         sigma(6) = sigma(6) + (1+3*h(r)+h(s))/((1+2*h(r))*(1+h(r)+h(s)));
      end
   end
end

for r = 1:p
   for s = 1:p
      for t = 1:p
         if (s~=r) && (t~=r) && (t~=s)
            sigma(7) = sigma(7) + (1+2*h(r)+h(s)+h(t))/((1+h(r)+h(s))*(1+h(r)+h(t)));
         end
      end
   end
end

rho4 = 2/d*(6*sigma(4)+3*sigma(5)+4*sigma(6)+sigma(7));

T = rho4/8 - 1/24*(3*rho13sq + 2*rho23sq);

   