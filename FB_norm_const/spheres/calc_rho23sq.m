function out = calc_rho23sq(CHat,p,dcum,mu,V)

runSum = 0;

lambda = V\mu;
K2HatInv = inv(calc_K2(CHat,p,dcum,mu,V));
CHatInv = inv(CHat);

K3Hat = zeros(p,p,p);

for k1 = 1:p
   for k2 = k1:p
      for k3 = k2:p
         K3Hat(k1,k2,k3) = calc_KHat_n(3,CHatInv,p,dcum,[k1 k2 k3],lambda);
      end
   end
end

for i = 1:p
   for j = 1:p
      for k = 1:p
         for r = 1:p
            for s = 1:p
               for t = 1:p
                  inda = sort([i j k]);
                  indb = sort([r s t]);
                  runSum = runSum + ...
                     K3Hat(inda(1),inda(2),inda(3))*...
                     K3Hat(indb(1),indb(2),indb(3))*...
                     K2HatInv(i,r)*K2HatInv(j,s)*K2HatInv(k,t);
               end
            end
         end
      end
   end
end

out = runSum;