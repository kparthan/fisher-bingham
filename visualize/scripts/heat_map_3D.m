function [] = heat_map_3D(file_name)

M = load(file_name);

x = M(:,1);
y = M(:,2);
z = M(:,3);
density = M(:,4);



h=scatter3(x,y,z,'cdata',density);

view ([90 180]);

%set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
%set(gca,'XTickLabel','');
%set(gca,'YTickLabel','');
%set(gca,'ZTickLabel','');

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

%saveas(gcf,'../ecc_1.eps','epsc');

end
