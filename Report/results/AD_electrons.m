omega = [0.01 0.1 0.28 0.5 0.75 1];
N2 = [14.7 3.36 1.76 1.24 9.62e-1 8.17e-1];

N6 = [99.1 22.9 12.1 8.56 6.7 5.57];

subplot(2,2,1)
plot(omega,N2,'o-');
xlabel('omega');
ylabel('Average distance');


subplot(2,2,3)
plot(omega,N6,'o-');
xlabel('omega');
ylabel('Average distance')

subplot(1,2,2)
plot(omega,N6/N6(1),'ro-',omega,N2/N2(1),'bo--');
xlabel('Omega');
ylabel('Average distance rel. to AD(0.01)')
legend('6 electrons','2 electrons')