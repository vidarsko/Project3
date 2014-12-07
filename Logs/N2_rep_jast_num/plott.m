
clear
E = csvread('energies.csv');
V = csvread('variances.csv');
alpha = csvread('alpha.csv');
beta = csvread('beta.csv');
subplot(2,1,1)
surf(beta,alpha,log10(V))   ;
xlabel('beta')
ylabel('alpha')
zlabel('log10(Variance)')

subplot(2,1,2)
surf(beta,alpha,E);
xlabel('beta')
ylabel('alpha')
zlabel('Energy')