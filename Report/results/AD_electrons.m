omega = [0.01 0.1 0.28 0.5 0.75 1];
N2 = [14.7 3.36 1.76 1.24 9.62e-1 8.17e-1];

N6 = [39.6 9.16 4.83 3.42 2.68 2.23];

subplot(2,2,1)
plot(omega,N2,'bo--');
xlabel('omega');
ylabel('AD');


subplot(2,2,3)
plot(omega,N6,'ro-');
xlabel('omega');
ylabel('AD')

subplot(1,2,2)
plot(omega,N6/N6(1),'ro-',omega,N2/N2(1),'bo--');
xlabel('Omega');
ylabel('$\frac{AD}{AD(0.01}$')
legend('6 electrons','2 electrons')