function [] = plot_prior_density()

clc;

k_inc = 0.025;
e_inc = 0.0005;
kappa = [0:k_inc:5-k_inc];
ecc = [0:e_inc:1-e_inc];

num_k = size(kappa,2)
num_e = size(ecc,2)
M = zeros(num_k,num_e);

fig = figure();
hold on;

for i = 1:num_k
  k = kappa(i);
  ksq_term = (1+k*k)^2;
  hk = (8/pi) * k / ksq_term;
  for j = 1:num_e
    e = ecc(j);
    %he = 6 * e * (1-e); % B(2,2)
    %he = 30 * e * (1-e)^4; % B(2,5)
    he = 30 * e^4 * (1-e); % B(5,2)
    %he = 110 * e * (1-e)^9; % B(2,10)
    %he = 110 * e^9 * (1-e); % B(10,2)
    M(i,j) = hk * he;
  end
end

mesh(M);
xlabel('X');
ylabel('Y');
zlabel('Z');
