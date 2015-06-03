function [] = heat_map_2D(file_name)

  addpath('export_fig');

  %[theta phi] = meshgrid(0:1:359.9,0:1:179.9);
  theta = 0:1:179.9;
  phi = 0:1:359.9;

  fig = figure();
  hold on;
  set(gcf, 'Color', 'w');
  xlabel('Longitude\phi','fontsize',12);
  ylabel('Co-latitude\theta','fontsize',13);
  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  set(xlabh,'interpreter','tex');
  set(ylabh,'interpreter','tex');
  set(gca,'Ylim',[15 115]);
  set(gca,'Xlim',[0 360]);
  set(gca,'xtick',[0:60:360],'fontsize',10);
  set(gca,'ytick',[15:20:180],'fontsize',10);
  %view ([0 90]);
  view ([-61 28]);

  M = load(file_name);
  %mesh(M);
  %surf(M);
  %shading interp;
  
  surf(M,'EdgeColor','none','LineStyle','none');
  %bar3(M);
  %colormap(jet);

%  xlabel('\phi');
%  ylabel('\theta');
%  zlabel('Z');
%
%  set(gca,'XLim',[0 360]);
%  set(gca,'YLim',[15 115]);

  %outfile = 'b_vmf_37_density';
  outfile = 'b_empirical';
  output_fig = strcat('../figs/protein_modelling/',outfile,'.fig');
  output_pdf = strcat('../figs/protein_modelling/',outfile,'.pdf');
  output_eps = strcat('../figs/protein_modelling/',outfile,'.eps');

  %saveas(gcf,output_fig);
  %export_fig(output_eps,'-eps');
  export_fig(output_pdf,'-pdf');

end
