function C = calc_C(theta,p,dcum,V)

C = inv(V);

for i = 1:p
   C = C - 2*theta(i)*calc_J(i,dcum);
end