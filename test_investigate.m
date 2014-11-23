E = csvread('energies.csv');
V = csvread('variances.csv');
alpha = csvread('alpha.csv');
beta = csvread('beta.csv');

subplot(2,2,1)
surf(beta,alpha,log10(V));
xlabel('beta')
ylabel('alpha')
zlabel('log(Variance)')

subplot(2,2,2)
surf(beta,alpha,V);
xlabel('beta')
ylabel('alpha')
zlabel('Variance')

subplot(2,2,3)
surf(beta,alpha,E);
xlabel('beta')
ylabel('alpha')

zlabel('Energy')


subplot(2,2,4)
surf(beta,alpha,log10(E));
xlabel('beta')
ylabel('alpha')
zlabel('log(Energy)')