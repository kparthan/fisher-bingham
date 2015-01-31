function [] = visualize()
  clc;
  D = 2;

  M = load('mixture_density.dat');
  fig = figure();

  hold on;
  scatter3(M(:,1),M(:,2),M(:,3),2,'cdata',M(:,3));
  set(gca,'XTick',[-4:2:4]);
  set(gca,'FontSize',15);
  xlabel('X');
  ylabel('Y');
  zlabel('Z');

  black = [0 0 0];
  red = [1 0 0];
  green = [0 0.5 0];
  blue = [0 0 1];

  % ... plotting the full/parent mixture
%  K = 4; mix_file = 'iter_3_c1_children_after_EM';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  split = 1;
%  for k=1:K
%    elipsnorm(mus(:,k),covs(:,:,k),2,0,black);
%  end

  % ... splitting (plotting the parent mixture)
%  K = 3; mix_file = 'iter_3_parent';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  split = 1;
%  for k=1:K
%    if (k ~= split)
%      elipsnorm(mus(:,k),covs(:,:,k),2,1,black);
%    end
%  end

  % ... splitting (plotting the initialization)
%  sp = [1.445, 2.026]; ep = [-1.349, 2.014];
%  draw_line(sp,ep);
%  status = 'iter_3_split_c1';
%  fig_file = strcat('./plots/',status,'.fig');
%  eps_file = strcat('./plots/',status,'.eps');

  % ... splitting (plotting the two children before/after EM adjustment)
%  K = 4; mix_file = 'iter_3_c1_children_before_EM';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  split1 = 1; split2 = 2;
%  elipsnorm(mus(:,split1),covs(:,:,split1),2,1,red);
%  elipsnorm(mus(:,split2),covs(:,:,split2),2,1,red);

  % ... deleting (plotting the parent mixture)
%  K = 3; mix_file = 'iter_3_parent';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  delete = 1;
%  for k=1:K
%    if (k == delete)
%      elipsnorm(mus(:,k),covs(:,:,k),2,0,green);
%    else
%      elipsnorm(mus(:,k),covs(:,:,k),2,0,black);
%    end
%  end
%  status = 'iter_3_delete_c1';

  % ... deleting (plotting the remaining components before/after EM adjustment)
%  K = 2; mix_file = 'iter_3_delete_c1_after_EM';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  for k=1:K
%    %elipsnorm(mus(:,k),covs(:,:,k),2,1,black); % before
%    elipsnorm(mus(:,k),covs(:,:,k),2,0,black);  % after
%  end

  % ... merging (plotting the parent mixture)
%  K = 3; mix_file = 'iter_3_parent';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  merge1 = 1; merge2 = 2;
%  for k=1:K
%    if (k == merge1 || k == merge2)
%      elipsnorm(mus(:,k),covs(:,:,k),2,0,blue);
%    else
%      elipsnorm(mus(:,k),covs(:,:,k),2,0,black);
%    end
%  end
%  status = 'iter_3_merge_c1_c2';

  % ... merging (plotting the merged components before EM adjustment)
  K = 3; mix_file = 'iter_3_parent';
  mix_path = strcat('./mixtures/',mix_file);
  [mus,covs] = parse_mixture_file(K,D,mix_path);
  merge1 = 1; merge2 = 2;
  for k=1:K
    if (k ~= merge1 && k ~= merge2)
      elipsnorm(mus(:,k),covs(:,:,k),2,1,black);
    end
  end
  K = 1; mix_file = 'iter_3_merge_c1_c2_before_EM';
  mix_path = strcat('./mixtures/',mix_file);
  [mus,covs] = parse_mixture_file(K,D,mix_path);
  elipsnorm(mus(:,1),covs(:,:,1),2,1,blue); % before

  % ... merging (plotting the merged components after EM adjustment)
%  K = 2; mix_file = 'iter_3_merge_c1_c2_after_EM';
%  mix_path = strcat('./mixtures/',mix_file);
%  [mus,covs] = parse_mixture_file(K,D,mix_path);
%  for k=1:K
%    elipsnorm(mus(:,k),covs(:,:,k),2,0,black); % before
%  end

  status = mix_file;
  fig_file = strcat('./plots/',status,'.fig');
  eps_file = strcat('./plots/',status,'.eps');

  savefig(fig,fig_file);
  saveas(fig,eps_file,'eps2c');
  %print(fig,'-deps2c',eps_file,'-r500');

end

% draw a line between sp and ep
function [] = draw_line(sp,ep)

  %set(gcf,'visible','off');
  x = [sp(1) ep(1)];
  y = [sp(2) ep(2)];
  plot(x,y,'--k','LineWidth',3);
  %plot(x,y,'o','Color',[0.3 0.3 0.3],'MarkerSize',10,'MarkerFaceColor',[0.3,0.3,0.3]);
  plot(x,y,'o','Color',[0.0 0.0 0.0],'MarkerSize',10,'MarkerFaceColor',[0.0,0.0,0.0]);

end

function y = elipsnorm(m,covar,level,dashed,col);
  if nargin<4 
     dashed=0;
  end
  [uu,ei,vv]=svd(covar);
  a = sqrt(ei(1,1)*level*level);
  b = sqrt(ei(2,2)*level*level);
  theta = [0:0.01:2*pi];
  xx = a*cos(theta);
  yy = b*sin(theta);
  cord = [xx' yy']';
  cord = uu*cord;
  if dashed==1
     dashline(cord(1,:)+m(1),cord(2,:)+m(2),1.5,3,1.5,3,'Color',col,'LineWidth',3);
  else
     plot(cord(1,:)+m(1),cord(2,:)+m(2),'Color',col,'LineWidth',3)
  end
end
