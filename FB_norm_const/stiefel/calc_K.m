function K = calc_K(theta,p,d,mu,V)

C = calc_C(theta,p,d,V);

K = -.5*log(det(C)) - .5*log(det(V)) + ...
  .5*(mu'/(V*C*V)*mu) -.5*(mu'/(V)*mu);