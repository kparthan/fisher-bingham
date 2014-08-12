function out = randVonMises(mean, k, n)

% RANDVONMISES  Draw a random sample from the von Mises distribution on the
% unit circle.  
% 
% Parameterised in terms of angle theta, the code returns n values of 
% theta from the distribution with density const*exp(k*cos(theta-mean)).
%
% This function is a direct transcoding of the R function rvm in the 
% package CircStats, itself an implementation of the algorithm in [1].
%
% Note that the function randVonMisesFisher.m covers the more general case
% of the von Mises(-Fisher) distribution in abritrary dimension, but is
% typically slower for the unit-circle case covered by this function.
%
% Reference:
% [1] Best & Fisher (1979) Efficient Simulation of the von Mises 
%     Distribution, JRSSC.

% Simon Preston (http://www.maths.nott.ac.uk/~spp), 2012

if nargin<3, n = 1; end

out = NaN(n,1);
a = 1 + (1 + 4 * (k^2))^0.5;
b = (a - (2 * a)^0.5)/(2 * k);
r = (1 + b^2)/(2 * b);

obs = 1;
while obs <= n
   U1 = rand;
   z = cos(pi * U1);
   f = (1 + r * z)/(r + z);
   c = k * (r - f);
   U2 = rand;
   if (c * (2 - c) - U2 > 0) 
      U3 = rand;
      out(obs) = sign(U3 - 0.5) * acos(f) + mean;
      out(obs) = mod(out(obs),(2 * pi));
      obs = obs + 1;
   else
      if (log(c/U2) + 1 - c >= 0)
         U3 = rand;
         out(obs) = sign(U3 - 0.5) * acos(f) + mean;
         out(obs) = mod(out(obs),(2 * pi));
         obs = obs + 1;
      end
   end
end      
