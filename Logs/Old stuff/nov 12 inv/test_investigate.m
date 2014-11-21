E = csvread('test_energies.csv');
V = csvread('test_variances.csv');
alpha = csvread('test_alpha.csv');
beta = csvread('test_beta.csv');

fig1 = surf(beta,alpha,log(V));
xlabel('beta')
ylabel('alpha')
zlabel('log(Variance)')
fig2u3d