function [] = visualize3D(K)

  % draw a unit sphere
  n = 25;
  r = ones(n, n); % radius is 1 
  [thetas, phis] = meshgrid(linspace(0, pi, n), linspace(0, 2*pi, n));

  X = zeros(n,n);
  Y = zeros(n,n);
  Z = zeros(n,n);
  for i=1:n
    for j=1:n
      theta = thetas(i,j);
      phi = phis(i,j);
      [x,y,z] = spherical2cartesian(theta,phi);
      X(i,j) = x;
      Y(i,j) = y;
      Z(i,j) = z;
    end
  end
%  [x,y,z]= sph2cart(th, phi, r);
  surface(X,Y,Z,'FaceColor','none');
  hold on;

  %c(1,:) = [1 0.5 0];

  % plot the sampled data
  for k = 1:K
     data_file = strcat('../sampled_data/comp',num2str(k),'.dat');
     M = load(data_file);
     x = M(:,1);
     y = M(:,2);
     z = M(:,3);
     colors = rand(1,3);
     %plot3(x,y,z,'.','Color',c(k,:));
     plot3(x,y,z,'.','Color',colors);
  end  
  xlabel('X');
  ylabel('Y');
  zlabel('Z');

  % create legend
  %N = [1:K];
  %legend_cell = cellstr('unit sphere');
  %legend_cell = [legend_cell ; cellstr(num2str(N','%d'))];
  %legend(legend_cell);

end
