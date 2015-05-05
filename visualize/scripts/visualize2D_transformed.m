function [] = visualize2D_transformed(K)

hold on;

% plot the sampled data
for k = 1:K
   %k
   data_file = strcat('../sampled_data/transformed_comp',num2str(k),'.dat');
   M = load(data_file);
   x = M(:,1);
   y = M(:,2);
   colors = rand(1,3);
   plot(x,y,'.','Color',colors);
end  

% create legend
N = [1:K];
legend_cell = [cellstr(num2str(N','%d'))];
%legend(legend_cell);

end
