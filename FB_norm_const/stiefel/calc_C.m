function C = calc_C(theta_m,p,d,V)

% Note the following dimensions:
%     dim(theta_m) = [p p];
%     dim(V) = [d*p d*p];
%     dim(C) = [d*p d*p];

C = inv(V);

for i = 1:p
   for j = i:p
      Jij = calc_J(i,j,p,d);
      C = C - theta_m(i,j)*(Jij + Jij');
   end
end