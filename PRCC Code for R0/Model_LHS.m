tic
clearvars;
close all;
clc;

R_0 = 1.726985;

runs=100000;
           
Parameter_settings_LHS;
beta_LHS=LHS_Call(0.9*0.9,beta,0.9*1.1,0,runs,'unif');
mu_LHS=LHS_Call(0.0092*0.9,mu,0.0092*1.1,0,runs,'unif');
gamma_LHS=LHS_Call(0.3*0.9,gamma,0.3*1.1,0,runs,'unif');
alpha_LHS=LHS_Call(0.7384*0.9,alpha,0.7384*1.1,0,runs,'unif');
psi_LHS=LHS_Call(0.384*0.9,psi,0.384*1.1,0,runs,'unif');
rho_LHS=LHS_Call(0.5*0.9,rho,0.5*1.1,0,runs,'unif');
phi_LHS=LHS_Call(0.6*0.9,phi,0.6*1.1,0,runs,'unif');
theta_LHS=LHS_Call(0.5*0.9,theta,0.5*1.1,0,runs,'unif');
sigma_LHS=LHS_Call(0.4*0.9,sigma,0.4*1.1,0,runs,'unif');
omega_LHS=LHS_Call(0.1*0.9,omega,0.1*1.1,0,runs,'unif');

LHSmatrix=[beta_LHS,...
    mu_LHS,...
    gamma_LHS,...
    alpha_LHS,...
    psi_LHS,...
    rho_LHS,...
    phi_LHS,...
    theta_LHS,...
    sigma_LHS,...
    omega_LHS];


for x=1:runs 
    x;
    LHSmatrix (x,:);
    R0 = RRR(LHSmatrix,x);
    R0_lhs(:,x)= [R0];
end

E_R = mean(R0_lhs);
disp(['E[R] = ', num2str(E_R)]);
disp(['R_0 = ', num2str(R_0)]);

figure(2);
histogram(R0_lhs, 50, 'Normalization', 'pdf');
hold on;
xlabel('Distibution of $\mathcal{R}_0$');
ylabel('Probability Density');
xline(E_R, 'Color', 'r', 'LineWidth', 2)
xline(R_0, 'Color', 'k', 'LineWidth', 2)
xline(prctile(R0_lhs, 5), 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2)%, 'Label', 'R_0(5th)');
xline(prctile(R0_lhs, 95), 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2)%, 'Label', 'R_0(95th)')
set(gca, 'LineWidth', 2);

%save Model_LHS.mat;

  [prcc sign sign_label]=PRCC(LHSmatrix,R0_lhs,1,PRCC_var,0.05);

toc