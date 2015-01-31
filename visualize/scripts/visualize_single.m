function [] = visualize_single(data_file)

% draw a unit sphere
n = 25;
r = ones(n, n); % radius is 1 
[th, phi] = meshgrid(linspace(0, pi, n), linspace(0, 2*pi, n));
[x,y,z]= sph2cart(th, phi, r);
surface(x,y,z,'FaceColor','none');
hold on;

% plot the sampled data
M = load(data_file);
x = M(:,1);
y = M(:,2);
z = M(:,3);
%colors = rand(1,3);
colors = [0 0 1];
plot3(x,y,z,'.','Color',colors);

xlabel('X');
ylabel('Y');
zlabel('Z');

% create legend
%N = [1:K];
%legend_cell = cellstr('unit sphere');
%legend_cell = [legend_cell ; cellstr(num2str(N','comp%-d'))];
%legend(legend_cell);

end
