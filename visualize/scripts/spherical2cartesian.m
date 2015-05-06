% unit sphere
function [x,y,z] = spherical2cartesian(theta,phi)

  x = cos(theta);
  y = sin(theta) * cos(phi);
  z = sin(theta) * sin(phi);

end
