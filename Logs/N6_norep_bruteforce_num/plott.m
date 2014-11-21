alpha = csvread('alpha.csv');
energies = csvread('energies.csv');
variances = csvread('variances.csv');
start = 1;

subplot(1,2,1);
plot(alpha(start:end),energies(start:end));
xlabel('alpha');
ylabel('Energies[a.u]');

subplot(1,2,2);
plot(alpha(start:end),log10(variances(start:end)));
xlabel('alpha');
ylabel('log10(variance)');