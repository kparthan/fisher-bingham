function [] = plot_posterior_prior2d()

  addpath('export_fig');

  %file_name = 'prior2d_posterior1';
  %file_name = 'prior2d_posterior2';
  %file_name = 'prior2d_posterior3';
  %file_name = 'prior2d_posterior2_1';
  %file_name = 'prior2d_posterior3_1';

  data_file = strcat('../sampled_data/',file_name,'.dat');
  output_fig = strcat('../figs/',file_name,'.fig');
  output_eps = strcat('../figs/',file_name,'.eps');
  output = strcat('../figs/',file_name);
  output_pdf = strcat('../figs/',file_name,'.pdf');
  output_jpg = strcat('../figs/',file_name,'.jpg');

  M = load(data_file);

  fig = figure();
  hold on;
  scatter3(M(:,1),M(:,2),M(:,3),'cdata',M(:,3));
  grid on;
  set(gcf, 'Color', 'w');

%  xlabel('\kappa','fontweight','bold','fontsize',25);
%  ylabel('\beta','fontweight','bold','fontsize',25);
%  set(gca,'xtick',[0:5:30],'fontsize',18);
%  set(gca,'ytick',[0:3:15],'fontsize',18);

% xlabel('\kappa','fontweight','bold','fontsize',25);
% ylabel('eccentricity (e)','fontsize',25);
% set(gca,'xtick',[0:5:30],'fontsize',18);
% set(gca,'ytick',[0:0.2:1],'fontsize',18);

%  xlabel('z_4','fontsize',25);
%  ylabel('z_5','fontsize',25);
%  set(gca,'xtick',[0:0.2:1],'fontsize',18);
%  set(gca,'ytick',[0:0.2:1],'fontsize',18);

%  saveas(gcf,output_fig);
%  saveas(fig,output_eps,'eps2c')

%  print2eps(output_eps);
%  eps2pdf(output_eps,output_pdf,1);

  %export_fig(output_pdf,'-pdf');
%  export_fig(output_jpg,'-jpg');

end
