function out = volumeOfStiefelManifold(d,p)

% The Stiefel manifold is the set {X : X'*X = eye(p)} where X is d-by-p.
%
% The implemented expression comes from the note "Volume of the Stiefel 
% Manifold" by Jason D. M. Rennie, 2006.

out = 1;

for k = (d-p+1):d
   out = out*2*pi^(k/2)/gamma(k/2);
end
