function [] = visualize_mixture_contours(K)

  bins_folder = '../sampled_data/bins_kent/';
%  bins_folder = '../sampled_data/bins_vmf/';
 
   % plot the entire mixture density
  data_file = strcat(bins_folder,'mixture_density.dat');
  M = load(data_file);
  density = M(:,4);
  n = size(M,1);
  angles = zeros(n,2);
  for i=1:n
    x = M(i,1);
    y = M(i,2);
    z = M(i,3);
    [phi, theta] = cartesian2spherical(x,y,z);
    angles(i,1) = phi;
    angles(i,2) = theta;
  end
  scatter3(angles(:,1),angles(:,2),density,2,'cdata',density);

  hold on;

  % plot the contours 
%  theta = 0:1:179.9;
%  phi = 0:1:359.9;
%  for k = 1:K
%    %k
%    data_file = strcat(bins_folder,'comp',num2str(k),'_prob_bins2D.dat');
%    prob_bins = load(data_file);
%    colors = rand(1,3);
%    %contour(phi,theta,prob_bins,[0 .1 .5 1],'LineWidth',2);
%    [C,h] = contour(phi,theta,prob_bins,1,'LineWidth',2);
%    %clabel(C,h);
%  end  

  % create legend
  %N = [1:K];
  %legend_cell = [cellstr(num2str(N','%d'))];
  %legend(legend_cell);

end
