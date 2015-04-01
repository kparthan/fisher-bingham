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

h=scatter3(x,y,z,'cdata',density);
%h=scatter3(x,y,z,2,'cdata',density);

%view ([180 90]);
%view ([0 90]);

%xlabel('X2');
%ylabel('X3');
%zlabel('Z');

%set(gca,'xtick',[-1,-0.5,0,0.5,1]);
%set(gca,'ytick',[-1,-0.5,0,0.5,1]);
%set(gca,'ztick',[]);

%saveas(gcf,'../figs/k10_e1_heatmap.fig');
%saveas(gcf,'../figs/k10_e1_heatmap.jpg');
%saveas(gcf,'../figs/k10_e1_scatter.fig');
%saveas(gcf,'../figs/k10_e1_scatter.jpg');

end
