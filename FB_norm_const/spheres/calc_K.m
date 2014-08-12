function K = calc_K(C,p,dcum,mu,V)

K = -.5*log(det(C)) - .5*log(det(V)) + ...
  .5*(mu'/(V*C*V)*mu) -.5*(mu'/(V)*mu);
