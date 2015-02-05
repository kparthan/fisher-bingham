function [] = visualize_single_transformed(data_file)

hold on;

% plot the sampled data
M = load(data_file);
x = M(:,1);
y = M(:,2);
%colors = rand(1,3);
colors = [0 0 1];
plot(x,y,'.','Color',colors);

xlabel('X');
ylabel('Y');

% create legend
%N = [1:K];
%legend_cell = cellstr('unit sphere');
%legend_cell = [legend_cell ; cellstr(num2str(N','comp%-d'))];
%legend(legend_cell);

end
