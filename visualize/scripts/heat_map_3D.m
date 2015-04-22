function [] = heat_map_3D(file_name)

  % draw a unit sphere
  %n = 25;
  %r = ones(n, n); % radius is 1 
  %[th, phi] = meshgrid(linspace(0, pi, n), linspace(0, 2*pi, n));
  %[x,y,z]= sph2cart(th, phi, r);
  %surface(x,y,z,'FaceColor','none');
  %hold on;

  % plot data
  M = load(file_name);
  x = M(:,1);
  y = M(:,2);
  z = M(:,3);
  density = M(:,4);

  h=scatter3(x,y,z,'cdata',density);

end
