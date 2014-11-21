%test compare
figure(2)
alpha = csvread('jas_rep_alpha.csv');
beta = csvread('jas_rep_beta.csv');
energies= csvread('jas_rep_energies.csv');
energies(2,183)= 3.047;
surf(alpha,beta,energies);
xlabel('alpha')
ylabel('beta')