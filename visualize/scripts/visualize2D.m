function [] = visualize2D(K)

hold on;

% plot the sampled data
for k = 1:K
   k
   data_file = strcat('comp',num2str(k),'.dat');
   M = load(data_file);
   x = M(:,1);
   y = M(:,2);
   z = M(:,3);
   [phi,theta,r] = cart2sph(x,y,z);
   rows = size(M,1);
   for i=1:rows
       if phi(i) < 0
           phi(i) = 2*pi + phi(i);
       end    
       theta(i) = (pi/2) - theta(i);
   end
   angles = [phi theta] .* (180/pi);
   colors = rand(1,3);
   plot(angles(:,1),angles(:,2),'.','Color',colors);
end  

% create legend
N = [1:K];
legend_cell = [cellstr(num2str(N','comp%-d'))];
%legend(legend_cell);

end
