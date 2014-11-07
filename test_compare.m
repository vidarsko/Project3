%test compare
alpha = csvread('test_alpha.csv');
beta = csvread('test_beta.csv');
energy_diff = csvread('energy_diff.csv');
surf(beta,alpha,energy_diff);