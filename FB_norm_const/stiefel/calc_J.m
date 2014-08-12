function J = calc_J(i,j,p,d)

J = zeros(d*p);

J(1+(i-1)*d:d*i, 1+(j-1)*d:d*j) = eye(d);