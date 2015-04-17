function [] = plot_posterior_prior3d()

%  M = load('../sampled_data/prior3d_posterior1.dat');
%  fig = figure();
%  hold on;
%
%  scatter3(M(:,1),M(:,2),M(:,3),2,'cdata',M(:,3));
%
%  xlabel('\kappa','fontweight','bold','fontsize',20);
%  ylabel('\beta','fontweight','bold','fontsize',20);
%  zlabel(['posterior',10,'density'],'fontweight','bold');
%
%  xlabh = get(gca,'XLabel');
%  set(xlabh,'Position',get(xlabh,'Position') + [0 -1 0])
%  ylabh = get(gca,'YLabel');
%  set(ylabh,'Position',get(ylabh,'Position') + [24 0 0])
%  zlabh = get(gca,'ZLabel');
%  set(zlabh,'Position',get(zlabh,'Position') + [-3 0 0.25e-6])
%
%  view ([42 75]);
%
%  grid on;
%  set(gca,'xtick',[0:4:20]);
%  set(gca,'ytick',[0:2:10]);
%  saveas(gcf,'../figs/prior3d_posterior1.fig');
%  saveas(gcf,'../figs/prior3d_posterior1.jpg');

  M = load('../sampled_data/prior3d_posterior2.dat');
  fig = figure();
  hold on;

  scatter3(M(:,1),M(:,2),M(:,3),2,'cdata',M(:,3));

  xlabel('\kappa','fontweight','bold','fontsize',20);
  ylabel('e','fontweight','bold','fontsize',16);
  zlabel(['posterior',10,'density'],'fontweight','bold');

  xlabh = get(gca,'XLabel');
  set(xlabh,'Position',get(xlabh,'Position') + [0 -0.1 0])
  ylabh = get(gca,'YLabel');
  set(ylabh,'Position',get(ylabh,'Position') + [24 0 0])
  zlabh = get(gca,'ZLabel');
  set(zlabh,'Position',get(zlabh,'Position') + [-3 0 1e-6])

  view ([42 75]);

  grid on;
  set(gca,'xtick',[0:4:20]);
  set(gca,'ytick',[0:0.2:1]);
  set(gca,'ztick',[0:2.5e-6:7.5e-6]);
  saveas(gcf,'../figs/prior3d_posterior2.fig');
  saveas(gcf,'../figs/prior3d_posterior2.jpg');

end
