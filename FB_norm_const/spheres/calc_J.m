function J = calc_J(i,dcum)

J = zeros(dcum(end));

J(1+dcum(i):dcum(i+1), 1+dcum(i):dcum(i+1)) = eye(dcum(i+1)-dcum(i));