function K2 = calc_K2(theta,p,d,mu,V)

leng = 1/2*p*(p+1);
K2 = zeros(leng,leng);
C = calc_C(theta,p,d,V);

for k1 = 1:leng
    for k2 = 1:leng
        [r1,s1] = get_r_s(k1);
        [r2,s2] = get_r_s(k2);
        Jr1s1 = calc_J(r1,s1,p,d);
        Jr2s2 = calc_J(r2,s2,p,d);
        K2(k1,k2) = ...
            1/2 * trace(C\(Jr2s2+Jr2s2')/C*(Jr1s1+Jr1s1')) + ...
            1/2 * mu'/(C*V)*(Jr2s2+Jr2s2')/C*(Jr1s1+Jr1s1')/(V*C)*mu + ...
            1/2 * mu'/(C*V)*(Jr1s1+Jr1s1')/C*(Jr2s2+Jr2s2')/(V*C)*mu;  
    end
end