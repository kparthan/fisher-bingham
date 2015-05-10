function [] = heat_map_2D(file_name)

  %[theta phi] = meshgrid(0:1:359.9,0:1:179.9);
  theta = 0:1:179.9;
  phi = 0:1:359.9;

  M = load(file_name);

  mesh(M);
  %bar3(M);
  %colormap(jet);

  xlabel('Longitude');
  ylabel('Co-latitude');
  zlabel('Z');

  set(gca,'XLim',[0 360]);
  %set(gca,'YLim',[0 125]);

end
