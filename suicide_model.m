%% ============================================================
%  Simulation of the suicide model and computation of R0
%% ============================================================

clear; clc; close all;

%% -----------------------------
% Parameters
%% -----------------------------
Lambda = 2000;
mu     = 0.0092;
beta   = 0.9000;
eta    = 1.0000;
gamma  = 0.300;
tau    = 0.3840;   % using psi in the table as tau
alpha  = 0.7384;
rho    = 0.5000;
theta  = 0.5;
sigma  = 0.400;
omega  = 0.100;
phi    = 0.600;

%% -----------------------------
% Initial conditions
% x = [S, M, I, Af, As, T, R]
%% -----------------------------
S0  = 99980;
M0  = 1;
I0  = 4;
Af0 = 0;
As0 = 2;
T0  = 1;
R0i = 1;

x0 = [S0; M0; I0; Af0; As0; T0; R0i];

%% -----------------------------
% Time span
%% -----------------------------
tspan = [0 50];   % years
tEval = linspace(tspan(1), tspan(2), 1000);

%% -----------------------------
% Compute R0
%% -----------------------------
R0 = (beta*gamma*((theta+mu) + (1-rho)*alpha)*(sigma+mu)*(omega+mu)) / ...
     ((alpha+tau+mu)*(theta+mu)*(gamma+mu)*(sigma+mu)*(omega+mu) ...
     - gamma*(1-phi)*omega*sigma*(tau*(theta+mu) + theta*(1-rho)*alpha));

fprintf('Basic reproduction number R0 = %.6f\n', R0);

%% -----------------------------
% Solve ODE system
%% -----------------------------
opts = odeset('NonNegative',1:7,'RelTol',1e-8,'AbsTol',1e-10);

[t,x] = ode15s(@(t,x) suicide(t,x,Lambda,mu,beta,eta,gamma,tau,alpha,rho,theta,sigma,omega,phi), ...
               tEval, x0, opts);

%% -----------------------------
% Display final solution
%% -----------------------------
fprintf('\nFinal state values at t = %.2f years:\n', t(end));
fprintf('S  = %.6f\n', x(end,1));
fprintf('M  = %.6f\n', x(end,2));
fprintf('I  = %.6f\n', x(end,3));
fprintf('Af = %.6f\n', x(end,4));
fprintf('As = %.6f\n', x(end,5));
fprintf('T  = %.6f\n', x(end,6));
fprintf('R  = %.6f\n', x(end,7));

%% -----------------------------
% Plot all state variables
%% -----------------------------
figure('Color','w');
%plot(t, x(:,1), 'b-',  'LineWidth', 2); 
hold on;
plot(t, x(:,2), 'r-',  'LineWidth', 2);
%plot(t, x(:,3), 'c-',  'LineWidth', 2);
%plot(t, x(:,4), 'k-',  'LineWidth', 2);
%plot(t, x(:,5), 'm-',  'LineWidth', 2);
%plot(t, x(:,6), 'g-',  'LineWidth', 2);
%plot(t, x(:,7), 'y-',  'LineWidth', 2);

xlabel('Time (years)');
ylabel('Population per 100,000');
title(sprintf('Suicide model dynamics (R_0 = %.4f)', R0));
%legend('M','I','A_s','Location','best');
box on;
grid on;

figure('Color','w');
hold on;
plot(t, x(:,3), 'c-',  'LineWidth', 2);
xlabel('Time (years)');
ylabel('Population per 100,000');
legend('I','Location','best');
box on;
grid on;

figure('Color','w');
hold on;
plot(t, x(:,5), 'k-',  'LineWidth', 2);

xlabel('Time (years)');
ylabel('Population per 100,000');
legend('A_s','Location','best');
box on;
grid on;

%% -----------------------------
% Plot living population only
%% -----------------------------
N = x(:,1) + x(:,2) + x(:,3) + x(:,5) + x(:,6) + x(:,7);

% figure('Color','w');
% plot(t, N, 'k-', 'LineWidth', 2.5);
% xlabel('Time (years)');
% ylabel('Living population per 100,000');
% title('Total living population');
% box on;
% grid on;

%% -----------------------------
% Plot cumulative suicide deaths
%% -----------------------------
% figure('Color','w');
% plot(t, x(:,4), 'r-', 'LineWidth', 2.5);
% xlabel('Time (years)');
% ylabel('Cumulative suicide deaths');
% title('Cumulative suicide-death class A_f(t)');
% box on;
% grid on;

%% ============================================================
% ODE function
%% ============================================================
function dx = suicide(~, x, Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi)

    % States
    S  = x(1);
    M  = x(2);
    I  = x(3);
    Af = x(4); %#ok<NASGU>
    As = x(5);
    T  = x(6);
    R  = x(7);

    % Living population
    N = S + M + I + As + T + R;
    if N <= 0
        N = 1e-12;
    end

    lambda = (I + As)/N;

    dx = zeros(7,1);

    dx(1) = Lambda - beta*lambda*S - mu*S + phi*omega*R;
    dx(2) = beta*lambda*S - (gamma + eta*I + mu)*M + (1-phi)*omega*R;
    dx(3) = (gamma + eta*I)*M - (alpha + tau + mu)*I;
    dx(4) = rho*alpha*I;
    dx(5) = (1-rho)*alpha*I - (theta + mu)*As;
    dx(6) = tau*I + theta*As - (sigma + mu)*T;
    dx(7) = sigma*T - (omega + mu)*R;
end