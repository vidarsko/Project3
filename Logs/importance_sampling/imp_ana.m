%importance analysis

dt = csvread('imp_ana_dt_meshgrid.csv');
MCS =csvread('imp_ana_MCS_meshgrid.csv');
BE = csvread('imp_ana_brute_energies.csv');
IE = csvread('imp_ana_imp_energies.csv');

% figure(1)
% hold off
% surf(log10(dt),MCS,BE);
% xlabel('log10(dt)');
% ylabel('MCS');
% zlabel('Brute energy');
% 
% figure(2)
% hold off
% surf(log10(dt),MCS,IE);
% xlabel('log10(dt)');
% ylabel('MCS');
% zlabel('Brute energy');


hold off
subplot(2,2,1)
i = 1;
plot(log10(dt(i,:)),BE(i,:),'bo-');
hold on
plot(log10(dt(i,:)),IE(i,:),'ro-');
xlabel('log10(dt)')
ylabel(strcat('Energy, MCS = 10^',num2str(log10(MCS(i,1)))));


subplot(2,2,2)
i = 2;
plot(log10(dt(i,:)),BE(i,:),'bo-');
hold on
plot(log10(dt(i,:)),IE(i,:),'ro-');
xlabel('log10(dt)')
ylabel(strcat('Energy, MCS = 10^',num2str(log10(MCS(i,1)))));

subplot(2,2,3)
i = 3;
plot(log10(dt(i,:)),BE(i,:),'bo-');
hold on
plot(log10(dt(i,:)),IE(i,:),'ro-');
xlabel('log10(dt)')
ylabel(strcat('Energy, MCS = 10^',num2str(log10(MCS(i,1)))));

subplot(2,2,4)
i = 4;
plot(log10(dt(i,:)),BE(i,:),'bo-');
hold on
plot(log10(dt(i,:)),IE(i,:),'ro-');
xlabel('log10(dt)')
ylabel(strcat('Energy, MCS = 10^',num2str(log10(MCS(i,1)))));
legend('Brute force energy','Importance sampling energy');


