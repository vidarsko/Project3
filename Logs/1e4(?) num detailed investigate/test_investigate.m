E = csvread('test_energies.csv');
V = csvread('test_variances.csv');
alpha = csvread('test_alpha.csv');
beta = csvread('test_beta.csv');
figure(1)
surf(beta,alpha,log(V));
xlabel('beta')
ylabel('alpha')
zlabel('log(Variance)')

figure(2)
surf(beta,alpha,E);
xlabel('beta')
ylabel('alpha')
zlabel('Energy')