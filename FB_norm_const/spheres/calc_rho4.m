function out = calc_rho4(CHat,p,dcum,mu,V)

runSum = 0;

lambda = V\mu;
K2HatInv = inv(calc_K2(CHat,p,dcum,mu,V));
CHatInv = inv(CHat);

K4Hat = zeros(p,p,p,p);

for k1 = 1:p
   for k2 = k1:p
      for k3 = k2:p
         for k4 = k3:p
            K4Hat(k1,k2,k3,k4) = ...
               calc_KHat_n(4,CHatInv,p,dcum,[k1 k2 k3 k4],lambda);
         end
      end
   end
end

for i = 1:p
   for j = 1:p
      for k = 1:p
         for l = 1:p
            ind = sort([i j k l]);
            runSum = runSum + ...
               K4Hat(ind(1),ind(2),ind(3),ind(4))*K2HatInv(i,j)*K2HatInv(k,l);
         end
      end
   end
end

out = runSum;