function [] = visualize3D(K)

% draw a unit sphere
n = 25;
r = ones(n, n); % radius is 1 
[th, phi] = meshgrid(linspace(0, pi, n), linspace(0, 2*pi, n));
[x,y,z]= sph2cart(th, phi, r);
surface(x,y,z,'FaceColor','none');
hold on;

c(1,:) = [1 0.5 0];
c(2,:) = [1 0 1];
c(3,:) = [0 1 1];
c(4,:) = [1 0 0];
c(5,:) = [0 1 0];
c(6,:) = [0 0 1];
c(7,:) = [0.5 0.5 0.5];
c(8,:) = [0.5 0.8 0.8];
c(9,:) = [0.25 0.25 0.25];

% plot the sampled data
for k = 1:K
   data_file = strcat('../comp',num2str(k),'.dat');
   M = load(data_file);
   x = M(:,1);
   y = M(:,2);
   z = M(:,3);
   colors = rand(1,3);
   plot3(x,y,z,'.','Color',c(k,:));
end  
xlabel('X');
ylabel('Y');
zlabel('Z');

% create legend
N = [1:K];
legend_cell = cellstr('unit sphere');
legend_cell = [legend_cell ; cellstr(num2str(N','%d'))];
%legend(legend_cell);

end