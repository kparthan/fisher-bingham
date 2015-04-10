function [] = plot_posterior_prior2d()

%  M = load('../sampled_data/prior2d_posterior1.dat');
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
%  saveas(gcf,'../figs/prior2d_posterior1.fig');
%  saveas(gcf,'../figs/prior2d_posterior1.jpg');

%  M = load('../sampled_data/prior2d_posterior2.dat');
%  fig = figure();
%  hold on;
%
%  scatter3(M(:,1),M(:,2),M(:,3),2,'cdata',M(:,3));
%
%  xlabel('\kappa','fontweight','bold','fontsize',20);
%  ylabel('e','fontweight','bold','fontsize',16);
%  zlabel(['posterior',10,'density'],'fontweight','bold');
%
%  xlabh = get(gca,'XLabel');
%  set(xlabh,'Position',get(xlabh,'Position') + [0 -0.1 0])
%  ylabh = get(gca,'YLabel');
%  set(ylabh,'Position',get(ylabh,'Position') + [24 0 0])
%  zlabh = get(gca,'ZLabel');
%  set(zlabh,'Position',get(zlabh,'Position') + [-3 0 1e-6])
%
%  view ([42 75]);
%
%  grid on;
%  set(gca,'xtick',[0:4:20]);
%  set(gca,'ytick',[0:0.2:1]);
%  saveas(gcf,'../figs/prior2d_posterior2.fig');
%  saveas(gcf,'../figs/prior2d_posterior2.jpg');

  M = load('../sampled_data/prior2d_posterior3.dat');
  fig = figure();
  hold on;

  scatter3(M(:,1),M(:,2),M(:,3),2,'cdata',M(:,3));

  xlabel('z_4','fontweight','bold','fontsize',16);
  ylabel('z_5','fontweight','bold','fontsize',16);
  zlabel(['posterior',10,'density'],'fontweight','bold');

  xlabh = get(gca,'XLabel');
  set(xlabh,'Position',get(xlabh,'Position') + [0 -0.1 0])
  ylabh = get(gca,'YLabel');
  set(ylabh,'Position',get(ylabh,'Position') + [1.2 0 0])
  zlabh = get(gca,'ZLabel');
  set(zlabh,'Position',get(zlabh,'Position') + [-0.2 0 0.01])

  view ([42 75]);

  grid on;
  %set(gca,'xtick',[0:4:20]);
  %set(gca,'ytick',[0:0.2:1]);
  saveas(gcf,'../figs/prior2d_posterior3.fig');
  saveas(gcf,'../figs/prior2d_posterior3.jpg');
end
