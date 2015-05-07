function [] = test()

  % plot the entire mixture density
  data_file = '../sampled_data/test.dat';
  M = load(data_file);

  x = M(:,1);
  y = M(:,2);
  z = M(:,3);

  min1 = min(z);
  max1 = max(z);
  range1 = max1 - min1;

  %hs = scatter3(x,y,z,1,'cdata',z);

  %set(gca,'XLim',[-30 30]);
  %set(gca,'YLim',[-30 30]);

  hold on;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');

  [xi, yi] = meshgrid(-5:0.1:15,-4:0.1:4);
  zi = griddata(x,y,z,xi,yi);
  %surf(xi,yi,zi);

  min2 = min(zi(:));
  max2 = max(zi(:));
  range2 = max2 - min2;

  factor = range1 / range2;
  zi2 = (zi-min2) * factor;
  zi3 = zi2 + min1;

  [C,h] = contour(xi,yi,zi3,'LineWidth',2,'LineColor','black');
  clabel(C,h);

end
