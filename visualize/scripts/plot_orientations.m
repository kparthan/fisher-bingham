function [] = plot_orientations()

  fig = figure();
  hold on;
  grid on;

  origin = [0,0,0];
  X2 = [1,0,0];
  X3 = [0,1,0];
  X1 = [0,0,1];
  X2neg = [-1,0,0];
  X3neg = [0,-1,0];
  X1neg = [0,0,-1];

  draw_line(origin,X1,'k','star4');
  draw_line(origin,X2,'k','star4');
  draw_line(origin,X3,'k','star4');
  draw_line(origin,X1neg,':k','star4');
  draw_line(origin,X2neg,':k','star4');
  draw_line(origin,X3neg,':k','star4');
  
  xlabel('X2');
  ylabel('X3');
  zlabel('X1');

  %savefig(fig,fig_file);
  %saveas(fig,eps_file,'eps2c');
end

% draw a line between sp and ep
function [] = draw_line(sp,ep,form1,form2)

  %set(gcf,'visible','off');
  x = [sp(1) ep(1)];
  y = [sp(2) ep(2)];
  z = [sp(3) ep(3)];
  plot3(x,y,z,form1,'LineWidth',3,'HeadStyle',form2);
  %plot3(x,y,z,'.','Color',[0 0 1]);
  %plot3(x,y,z,'o','Color',[0.3 0.3 0.3],'MarkerSize',10,'MarkerFaceColor',[0.3,0.3,0.3]);
  %plot(x,y,'o','Color',[0.0 0.0 0.0],'MarkerSize',10,'MarkerFaceColor',[0.0,0.0,0.0]);

end

