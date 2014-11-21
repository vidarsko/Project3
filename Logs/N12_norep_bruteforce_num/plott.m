alpha = csvread('N12alpha.csv');
energies = csvread('N12energies.csv');
variances = csvread('N12variances.csv');
start = 1;

subplot(1,2,1);
plot(alpha(start:end),energies(start:end));
xlabel('alpha');
ylabel('Energies[a.u]');

subplot(1,2,2);
plot(alpha(start:end),log10(variances(start:end)));
xlabel('alpha');
ylabel('log10(variance)');