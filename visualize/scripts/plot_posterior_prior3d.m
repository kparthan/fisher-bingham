function [] = plot_posterior_prior3d()

  addpath('export_fig');

  file_name = 'prior3d_posterior1';
  %file_name = 'prior3d_posterior2';
  %file_name = 'prior3d_posterior2_1';

  data_file = strcat('../sampled_data/',file_name,'.dat');
  output_fig = strcat('../figs/',file_name,'.fig');
  output_eps = strcat('../figs/',file_name,'.eps');
  output = strcat('../figs/',file_name);
  output_pdf = strcat('../figs/',file_name,'.pdf');

  M = load(data_file);

  fig = figure();
  hold on;
  scatter3(M(:,1),M(:,2),M(:,3),'cdata',M(:,3));
  grid on;
  set(gcf, 'Color', 'w');

  %view ([42 75]);
  view ([35 65]);

  xlabh = get(gca,'XLabel');
  ylabh = get(gca,'YLabel');
  zlabh = get(gca,'ZLabel');

  %set(xlabh,'interpreter','tex');
  %set(ylabh,'interpreter','tex');
  %zlabel(['posterior',10,'density'],'fontweight','bold');

  xlabel('\kappa','fontweight','bold','fontsize',25);
  ylabel('\beta','fontweight','bold','fontsize',25);
  set(xlabh,'Position',[15 -2 0]);
  set(ylabh,'Position',[35 7.5 0]);
  set(gca,'xtick',[0:5:30],'fontsize',12);
  set(gca,'ytick',[0:3:15],'fontsize',12);

% xlabel('\kappa','fontweight','bold','fontsize',25);
% ylabel('e','fontsize',25);
% set(xlabh,'Position',[15 -0.2 0]);
% set(ylabh,'Position',[35 0.5 0]);
% set(gca,'xtick',[0:5:30],'fontsize',12);
% set(gca,'ytick',[0:0.2:1],'fontsize',12);

  %saveas(gcf,output_fig);
%  saveas(fig,output_eps,'eps2c')

  %print2eps(output_eps);
  %eps2pdf(output_eps,output_pdf,1);

end
