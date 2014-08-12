function K1 = calc_K1(C,p,dcum,mu,V)

K1 = zeros(p,1);

for i = 1:p
   Ji = calc_J(i,dcum);
   K1(i) = trace(C\Ji) + mu'*((C*V)\Ji/(V*C))*mu;
end