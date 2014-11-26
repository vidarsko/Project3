%ds

parameters = csvread('N2_nojastrow.csv');
omega = parameters(:,1);
alpha = parameters(:,2);
beta = parameters(:,3);
energy = parameters(:,4);
variances = parameters(:,5);

figure(1)
plot(omega,alpha,'bo-')
hold on
plot(omega,beta,'ro-')

figure(2)
plot(omega,energy,'ro-');

figure(3)
plot(omega,variances,'bo-')
