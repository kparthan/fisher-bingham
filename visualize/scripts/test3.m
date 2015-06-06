function [] = test3(data_file)

M = load(data_file);

fig = figure();
%axis equal;
hold on;
set(gcf, 'Color', 'w');
scatter3(M(:,2),M(:,1),M(:,3),0.1,'cdata',M(:,3));
set(gca,'xlim',[0 360]);
%set(gca,'ylim',[0 180]);
set(gca,'Ylim',[15 115]);
set(gca,'xtick',[0:60:360]);
set(gca,'ytick',[15:20:180]);
xlabel('\phi');
ylabel('\theta');
view ([0 90]);

end
