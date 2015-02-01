function [] = stats(file_name)

  M = load(file_name);
  kappas = M(:,1);
  betas = M(:,2);
  ex = M(:,3);

  moment = M(:,4);
  ml = M(:,5);
  map = M(:,6);
  mml = M(:,7);

  figure();
  xlabel('ex');
  ylabel('kappas');
  zlabel('measure');

  hold on;

  X = 10:10:100;
  Y = 0.1:0.1:0.9;
  %meshgrid(10:10:100,0.1:0.1:0.9);

  Z_moment = vec2mat(moment,size(Y,2));
  surf(Y,X,Z_moment,'FaceColor',[1 0 0], 'FaceAlpha',0.5, 'EdgeAlpha', 0);

  Z_ml = vec2mat(ml,size(Y,2));
  surf(Y,X,Z_ml,'FaceColor',[0 1 0], 'FaceAlpha',0.5, 'EdgeAlpha', 0);

  Z_map = vec2mat(map,size(Y,2));
  surf(Y,X,Z_map,'FaceColor',[0 0 1], 'FaceAlpha',0.5, 'EdgeAlpha', 0);

  Z_mml = vec2mat(mml,size(Y,2));
  surf(Y,X,Z_mml,'FaceColor',[0 0 0], 'FaceAlpha',0.5, 'EdgeAlpha', 0);

  %colorbar;

end
