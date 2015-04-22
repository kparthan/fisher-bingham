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

  fig = figure();
  h=scatter3(x,y,z,'cdata',density);
  set(gcf, 'Color', 'w');

  %view ([180 90]);
  view ([0 90]);

  xlabel('X_2','fontsize',20);
  ylabel('X_3','fontsize',20);
  %zlabel('Z');

  set(gca,'xtick',[-1,-0.5,0,0.5,1],'fontsize',18);
  set(gca,'ytick',[-1,-0.5,0,0.5,1],'fontsize',18);
  %set(gca,'ztick',[]);

  file_name = 'k10_e9_heatmap';
  output_fig = strcat('../../figs/',file_name,'.fig');
  output_eps = strcat('../../figs/',file_name,'.eps');
  output_pdf = strcat('../../figs/',file_name,'.pdf');

  saveas(gcf,output_fig);

  print2eps(output_eps);
  eps2pdf(output_eps,output_pdf,1);

end
