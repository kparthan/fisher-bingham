function out = calc_KHat2_matrixFisher(phiHat_d,lambda_d,d)

% NB: phiHat_d is a vector of the diagonal elements of the (diagonal) phi

p = length(lambda_d);

leng = 1/2*p*(p+1);
out = zeros(leng,leng);

for i = 1:leng
    for j = 1:leng
        [r1,s1] = get_r_s(i);
        [r2,s2] = get_r_s(j);
        if r1==r2 && s1==s2
            if r2==s1
                out(i,j) =  2*d*phiHat_d(r1)^2 + 4*lambda_d(r1)^2*phiHat_d(r1)^3;
            else
                out(i,j) = d*phiHat_d(r1)*phiHat_d(s1) + ...
                    phiHat_d(r1)*phiHat_d(s1)*(lambda_d(r1)^2*phiHat_d(r1) +  ...
                    lambda_d(s1)^2*phiHat_d(s1));
            end
        end
    end
end
