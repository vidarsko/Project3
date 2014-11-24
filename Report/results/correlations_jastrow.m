%Correlations jastrow
hold off
omega = [0.01 0.1 0.28 0.5 0.75 1.0];
C = [0.282 0.16 0.105 0.078 0.064 0.048];

subplot(1,2,1);
plot(omega,C,'o-')
xlabel('omega');
ylabel('Correlation');

subplot(1,2,2);
plot(omega,log(C),'o-');
hold on
xlabel('omega')
ylabel('log(Correlation)')