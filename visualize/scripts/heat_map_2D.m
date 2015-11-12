function [] = heat_map_2D(file_name)

  addpath('export_fig');

  %[theta phi] = meshgrid(0:1:359.9,0:1:179.9);
  theta = 0:1:179.9;
  phi = 0:1:359.9;

  fig = figure();
  hold on;
  set(gcf, 'Color', 'w');
  xlabel('Longitude\phi','fontsize',12,'fontweight','bold');
  ylabel('Co-latitude\theta','fontsize',13,'fontweight','bold');
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  set(xlabh,'interpreter','tex');
  set(ylabh,'interpreter','tex');
  set(gca,'Xlim',[0 360]);
  set(gca,'xtick',[0:60:360],'fontsize',10);
  set(gca,'Ylim',[0 180]);
  set(gca,'ytick',[0:30:180],'fontsize',10);
  %view ([0 90]);
  view ([-27 18]);

  %set(xlabh,'Position',[150 -50 0]);
  %set(ylabh,'Position',[-15 90 0]);

  set(xlabh,'Rotation',5);
  set(ylabh,'Rotation',-25);

  M = load(file_name);
  mesh(M);
  %surf(M);
  %shading interp;
  
  %surf(M,'EdgeColor','none','LineStyle','none');
  %bar3(M);
  %colormap(jet);

%  xlabel('\phi');
%  ylabel('\theta');
%  zlabel('Z');
%
%  set(gca,'XLim',[0 360]);
%  set(gca,'YLim',[15 115]);

  %outfile = 'b_kent_23_2d';
  outfile = 'b_kent_23_3d';

  %outfile = 'b_empirical_2d';
  %outfile = 'b_empirical_3d';

  output_fig = strcat('../figs/arun/',outfile,'.fig');
  output_pdf = strcat('../figs/arun/',outfile,'.pdf');
  output_eps = strcat('../figs/arun/',outfile,'.eps');

  %saveas(gcf,output_fig);
  %export_fig(output_pdf,'-pdf');
  %print2eps(output_eps);

end
