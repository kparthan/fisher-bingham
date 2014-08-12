function K1 = calc_K1(theta_m,p,d,mu,V)

leng = 1/2*p*(p+1);
K1 = zeros(leng,1);
C = calc_C(theta_m,p,d,V);

for k = 1:leng
    [r,s] = get_r_s(k);
    Jrs = calc_J(r,s,p,d);
    K1(k) = 1/2 * trace(C\(Jrs + Jrs')) + ...
        1/2*mu'/(C*V)*(Jrs + Jrs')/(V*C)*mu;
end