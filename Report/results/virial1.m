%Virial theorem plots

omega = [0.01 0.1 0.28 0.5 0.75 1];

V2rep = [6.43e-2 3.5e-1 7.71e-1  1.22 1.67 2.12];
T2rep = [9.74e-3 9.09e-2  2.51e-1 4.44e-1  6.70e-1 8.81e-1];

V2norep = [9.94e-3   9.95e-2   2.78e-1 4.97e-1 7.46e-1 9.95e-1];
T2norep = [ 1.01e-2  1.01e-1   2.82e-1   5.03e-1 7.54e-1 1.01];

V6rep = [6.71e-1 3.24 6.68 10.1 13.5 16.6];
T6rep = [ 2.76e-2 3.30e-1 9.56e-1 1.74 2.68 3.62];

V6norep = [5.00e-2 4.99e-1  1.40 2.50 3.74 4.99];
T6norep = [5.00e-2 5.00e-1 1.40 2.50 3.76 5.00];


subplot(4,2,1)
plot(V2norep,T2norep,'bo-');
ylabel('Kinetic energy <T>')
xlabel('Potential energy <V>')

subplot(4,2,3)
plot(V2rep,T2rep,'ro-');
ylabel('Kinetic energy <T>')
xlabel('Potential energy <V>')


subplot(4,2,5)
plot(V6norep,T6norep,'mo-');
ylabel('Kinetic energy <T>')
xlabel('Potential energy <V>')

subplot(4,2,7)
plot(V6rep,T6rep,'ko-');
ylabel('Kinetic energy <T>')
xlabel('Potential energy <V>')

subplot(1,2,2)
plot(omega,T2norep./V2norep,'bo-', omega, T2rep./V2rep, 'ro-', omega, T6norep./V6norep,'mo-', omega, T6rep./V6rep,'ko-');
legend('N=2, no repulsion','N=2, with repulsion', 'N=6, no repulsion','N=6, with repulsion')
xlabel('omega')
ylabel('<T>/<V>')
