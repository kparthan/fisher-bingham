function [] = test2()

  x1 = -10:.02:10; x2 = -10:.02:10;
  [X1,X2] = meshgrid(x1,x2);
  F = mvnpdf([X1(:) X2(:)],[0 0],[9 0; 0 1]);
  F = reshape(F,length(x2),length(x1));
  F = F./sum(sum(F));

  %mesh(F);
  [FF,I] = sort(F(:),1,'descend');
  FFsum = cumsum(FF);
  % Load the values of FFsum into YY according to the original order of F
  % given by index I
  YY = ones(length(FF),1);
  YY(I) = FFsum;
  % Reshape YY to have the orignal form of F
  YYY = reshape(YY,size(F));
  % Plot the contours of the cumulative weight
%  figure
%  pcolor(x1,x2,F)
%  shading interp
  figure
  contour_levels = [[0.1:0.1:0.9] 0.86 0.98 0.99]; 
  contour_labels = [[0.1:0.2:0.9] 0.98 0.99];
  [C,h] = contour(x1,x2,YYY,contour_levels);
  clabel(C,h,contour_labels,'labelspacing',500)

  text(3,0,'+','Color','red');
  text(-3,0,'+','Color','red');
  text(6,0,'*','Color','red');
  text(-6,0,'*','Color','red');
  text(9,0,'x','Color','red');
  text(-9,0,'x','Color','red');

  %text(0,1,'*','Color','blue');
  %text(0,-1,'*','Color','blue');

end
