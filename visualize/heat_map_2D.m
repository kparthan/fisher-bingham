function [] = heat_map_2D(file_name)

%[theta phi] = meshgrid(0:0.1:359.9,0:0.1:179.9);
%M = load('heat_map.data');
M = load(file_name);
mesh(M);
%image([0 180], [0 360], M, 'CDataMapping', 'scaled');
%zlim([0,250]);