function [phi, theta] = cartesian2spherical(x,y,z)

  theta = acos(x);
  ratio = y/sin(theta);
  if (ratio > 1)
    ratio = 1;
  elseif ratio < -1
    ratio = -1;
  end
  angle = acos(ratio);
  phi = 0;
  if (y == 0 && z == 0)
    phi = 0;
  elseif (y == 0) 
    if (z > 0) 
      phi = angle;
    else
      phi = 2 * pi - angle;
    end
  elseif (z >= 0)
    phi = angle;
  elseif (z < 0)
    phi = 2 * pi - angle;
  end
  phi = phi * 180/pi;
  theta = theta * 180/pi;

end
