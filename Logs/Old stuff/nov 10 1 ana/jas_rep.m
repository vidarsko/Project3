%test compare
alpha = csvread('jas_rep_alpha.csv');
beta = csvread('jas_rep_beta.csv');
energies= csvread('jas_rep_energies.csv');
surf(beta,alpha,energies);
xlabel('beta')
ylabel('alpha')