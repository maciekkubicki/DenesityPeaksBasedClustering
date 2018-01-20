clear all;
clc
figure()
d = importdata('outputp.dat');
scatter(d(:,1), d(:,2));
xlabel('$\rho_i$','Interpreter','latex')
ylabel('$\delta_i$','Interpreter','latex')
