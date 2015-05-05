function [] = heat_map_2D(file_name)

%[theta phi] = meshgrid(0:1:359.9,0:1:179.9);
theta = 0:1:179.9;
phi = 0:1:359.9;
M = load(file_name);
%mesh(M);
meshc(phi,theta,M);
%contour(phi,theta,M,37);
xlabel('Longitude');
ylabel('Co-latitude');
zlabel('Z');

set(gca,'XLim',[0 360]);
set(gca,'YLim',[0 125]);
%set(gca,'ZLim',[-5 5]);
%image([0 180], [0 360], M, 'CDataMapping', 'scaled');
%zlim([0,250]);
