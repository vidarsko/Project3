%jastrow
alpha = csvread('jastrow_alpha.csv');
beta = csvread('jastrow_beta.csv');
energy =csvread('jastrow_energies.csv');

surf(beta,alpha,energy)
xlabel('beta')
ylabel('alpha')
