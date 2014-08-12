function K2 = calc_K2(C,p,dcum,mu,V)

K2 = zeros(p);

for i = 1:p
   for j = 1:p
      Ji = calc_J(i,dcum);
      Jj = calc_J(j,dcum);   
      K2(i,j) = 2*trace((C\Jj)*(C\Ji)) + 2*(mu')/(C*V)*(Ji/C*Jj + Jj/C*Ji)/(V*C)*mu;
   end
end