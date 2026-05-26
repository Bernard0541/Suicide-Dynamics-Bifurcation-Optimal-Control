function out = suicide_control
%% ============================================================
% Optimal control simulation for the suicide model
%
% FULL STABILIZED CODE
%
% Notes:
% 1. psi in the LaTeX is represented by tau in the MATLAB code.
% 2. This version lets you safely switch controls on/off.
% 3. It includes:
%    - smaller relaxation parameter
%    - shorter default horizon for stability
%    - finer time discretization
%    - NaN/Inf checks
%    - strategy flags
%    - separate state and control costs
%
% IMPORTANT:
% This code implements the adjoint system as currently written in your
% manuscript. If the adjoint equations are later re-derived to include
% the full N-dependence in the incidence term, then the code should be
% updated accordingly.
%% ============================================================

clearvars; clc; close all;

%% -----------------------------
% Strategy flags: switch controls on/off here
%% -----------------------------
use_xi1 = true;   % prevention control
use_xi2 = true;    % ideator-to-treatment enhancement
use_xi3 = true;    % non-fatal-attempter-to-treatment enhancement

%% -----------------------------
% Time settings
%% -----------------------------
tf     = 50;       % final time (years) -- reduced for stability
nSteps = 5000;     % finer grid
t      = linspace(0, tf, nSteps+1);
h      = tf / nSteps;

%% -----------------------------
% Parameters
%% -----------------------------
Lambda = 2000;
mu     = 0.0092;
beta   = 0.9000;
eta    = 1.0000;
gamma  = 0.3000;
alpha  = 0.7384;
rho    = 0.5000;
theta  = 0.5000;
sigma  = 0.4000;
omega  = 0.1000;
phi    = 0.6000;
tau    = 0.3840;   % psi in LaTeX = tau in code

%% -----------------------------
% Weight constants
%% -----------------------------
B1 = 120;
B2 = 150;
B3 = 100;

C1 = 100;
C2 = 100;
C3 = 100;

%% -----------------------------
% Initial conditions
%% -----------------------------
S0  = 99980;
M0  = 1;
I0  = 4;
Af0 = 0;
As0 = 2;
T0  = 1;
R0  = 1;

%% -----------------------------
% Allocate state variables
%% -----------------------------
S  = zeros(1,nSteps+1);
M  = zeros(1,nSteps+1);
I  = zeros(1,nSteps+1);
Af = zeros(1,nSteps+1);
As = zeros(1,nSteps+1);
Tt = zeros(1,nSteps+1);
R  = zeros(1,nSteps+1);

S(1)=S0; M(1)=M0; I(1)=I0; Af(1)=Af0; As(1)=As0; Tt(1)=T0; R(1)=R0;

%% -----------------------------
% Allocate adjoints
%% -----------------------------
lambda1 = zeros(1,nSteps+1);
lambda2 = zeros(1,nSteps+1);
lambda3 = zeros(1,nSteps+1);
lambda4 = zeros(1,nSteps+1);
lambda5 = zeros(1,nSteps+1);
lambda6 = zeros(1,nSteps+1);
lambda7 = zeros(1,nSteps+1);

%% -----------------------------
% Initial control guesses
%% -----------------------------
xi1 = 0.5*ones(1,nSteps+1);
xi2 = 0.5*ones(1,nSteps+1);
xi3 = 0.5*ones(1,nSteps+1);

if ~use_xi1, xi1(:) = 0; end
if ~use_xi2, xi2(:) = 0; end
if ~use_xi3, xi3(:) = 0; end

%% -----------------------------
% Iteration settings
%% -----------------------------
tol        = 1e-6;
maxIter    = 500;
relaxParam = 0.1;   % smaller relaxation for stability
err        = inf;
iter       = 0;

%% ============================================================
% Forward-backward sweep
%% ============================================================
while err > tol && iter < maxIter
    iter = iter + 1;

    % Store old values
    oldS = S; oldM = M; oldI = I; oldAf = Af; oldAs = As; oldTt = Tt; oldR = R;
    oldlambda1 = lambda1; oldlambda2 = lambda2; oldlambda3 = lambda3;
    oldlambda4 = lambda4; oldlambda5 = lambda5; oldlambda6 = lambda6; oldlambda7 = lambda7;
    oldxi1 = xi1; oldxi2 = xi2; oldxi3 = xi3;

    %% -----------------------------
    % Forward sweep with RK4
    %% -----------------------------
    for k = 1:nSteps
        xk = [S(k); M(k); I(k); Af(k); As(k); Tt(k); R(k)];
        uk = [xi1(k); xi2(k); xi3(k)];

        k1 = suicide_rhs_controlled(xk, uk, use_xi1, use_xi2, use_xi3, ...
            Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi);

        k2 = suicide_rhs_controlled(xk + 0.5*h*k1, uk, use_xi1, use_xi2, use_xi3, ...
            Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi);

        k3 = suicide_rhs_controlled(xk + 0.5*h*k2, uk, use_xi1, use_xi2, use_xi3, ...
            Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi);

        k4 = suicide_rhs_controlled(xk + h*k3, uk, use_xi1, use_xi2, use_xi3, ...
            Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi);

        xnext = xk + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        xnext = max(xnext, 0);

        S(k+1)  = xnext(1);
        M(k+1)  = xnext(2);
        I(k+1)  = xnext(3);
        Af(k+1) = xnext(4);
        As(k+1) = xnext(5);
        Tt(k+1) = xnext(6);
        R(k+1)  = xnext(7);
    end

    if any(~isfinite([S M I Af As Tt R]))
        error('State variables became non-finite during the forward sweep.');
    end

    %% -----------------------------
    % Backward sweep for adjoints
    %% -----------------------------
    lambda1(end)=0; lambda2(end)=0; lambda3(end)=0; lambda4(end)=0;
    lambda5(end)=0; lambda6(end)=0; lambda7(end)=0;

    for k = nSteps:-1:1
        Nk = S(k+1) + M(k+1) + I(k+1) + As(k+1) + Tt(k+1) + R(k+1);
        if Nk <= 0
            Nk = 1e-12;
        end

        xi1k = xi1(k+1); if ~use_xi1, xi1k = 0; end
        xi2k = xi2(k+1); if ~use_xi2, xi2k = 0; end
        xi3k = xi3(k+1); if ~use_xi3, xi3k = 0; end

        dl1 = beta*(lambda1(k+1)-lambda2(k+1))*(1-xi1k)*(I(k+1)+As(k+1))/Nk ...
              + mu*lambda1(k+1);

        dl2 = (lambda2(k+1)-lambda3(k+1))*(gamma + (1-xi1k)*eta*I(k+1)) ...
              + mu*lambda2(k+1) - B1;

        dl3 = alpha*(lambda3(k+1)-rho*lambda4(k+1)) ...
              - alpha*(1-rho)*lambda5(k+1) ...
              + (lambda3(k+1)-lambda6(k+1))*(1+xi2k)*tau ...
              + mu*lambda3(k+1) - B2 ...
              + (lambda2(k+1)-lambda3(k+1))*(1-xi1k)*eta*M(k+1) ...
              + beta*(lambda1(k+1)-lambda2(k+1))*(1-xi1k)*(I(k+1)+As(k+1))/Nk;

        dl4 = 0;

        dl5 = beta*(lambda1(k+1)-lambda2(k+1))*(1-xi1k)*S(k+1)/Nk ...
              + (lambda5(k+1)-lambda6(k+1))*(1+xi3k)*theta ...
              + mu*lambda5(k+1) - B3;

        dl6 = sigma*(lambda6(k+1)-lambda7(k+1)) + mu*lambda6(k+1);

        dl7 = (omega+mu)*lambda7(k+1) ...
              - (1-phi)*omega*lambda2(k+1) - phi*omega*lambda1(k+1);

        lambda1(k) = lambda1(k+1) - h*dl1;
        lambda2(k) = lambda2(k+1) - h*dl2;
        lambda3(k) = lambda3(k+1) - h*dl3;
        lambda4(k) = lambda4(k+1) - h*dl4;
        lambda5(k) = lambda5(k+1) - h*dl5;
        lambda6(k) = lambda6(k+1) - h*dl6;
        lambda7(k) = lambda7(k+1) - h*dl7;
    end

    if any(~isfinite([lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7]))
        error('Adjoint variables became non-finite during the backward sweep.');
    end

    %% -----------------------------
    % Control update
    %% -----------------------------
    Nall = S + M + I + As + Tt + R;
    Nall(Nall <= 0) = 1e-12;

    xi1_proj = (eta*(lambda3-lambda2).*I.*M.*Nall + beta*(lambda2-lambda1).*(I+As).*S) ./ (C1*Nall);
    xi2_proj = (tau*(lambda3-lambda6).*I) ./ C2;
    xi3_proj = (theta*(lambda5-lambda6).*As) ./ C3;

    xi1_proj = min(1, max(0, xi1_proj));
    xi2_proj = min(1, max(0, xi2_proj));
    xi3_proj = min(1, max(0, xi3_proj));

    if ~use_xi1
        xi1_proj(:) = 0;
    end
    if ~use_xi2
        xi2_proj(:) = 0;
    end
    if ~use_xi3
        xi3_proj(:) = 0;
    end

    xi1 = relaxParam*xi1_proj + (1-relaxParam)*oldxi1;
    xi2 = relaxParam*xi2_proj + (1-relaxParam)*oldxi2;
    xi3 = relaxParam*xi3_proj + (1-relaxParam)*oldxi3;

    if ~use_xi1, xi1(:) = 0; end
    if ~use_xi2, xi2(:) = 0; end
    if ~use_xi3, xi3(:) = 0; end

    %% -----------------------------
    % Convergence error
    %% -----------------------------
    err = max([
        norm(xi1-oldxi1, inf), norm(xi2-oldxi2, inf), norm(xi3-oldxi3, inf), ...
        norm(S-oldS, inf), norm(M-oldM, inf), norm(I-oldI, inf), ...
        norm(Af-oldAf, inf), norm(As-oldAs, inf), norm(Tt-oldTt, inf), norm(R-oldR, inf), ...
        norm(lambda1-oldlambda1, inf), norm(lambda2-oldlambda2, inf), ...
        norm(lambda3-oldlambda3, inf), norm(lambda4-oldlambda4, inf), ...
        norm(lambda5-oldlambda5, inf), norm(lambda6-oldlambda6, inf), ...
        norm(lambda7-oldlambda7, inf)
    ]);

    fprintf('Iteration %d, error = %.8e\n', iter, err);
end

if iter == maxIter && err > tol
    warning('Maximum iterations reached before convergence.');
end

%% -----------------------------
% Costs
%% -----------------------------
stateIntegrand   = B1*M + B2*I + B3*As;
controlIntegrand = 0.5*C1*xi1.^2 + 0.5*C2*xi2.^2 + 0.5*C3*xi3.^2;

stateCost   = trapz(t, stateIntegrand);
controlCost = trapz(t, controlIntegrand);
J = stateCost + controlCost;

%% -----------------------------
% Uncontrolled model on same grid
%% -----------------------------
Y0 = [S0; M0; I0; Af0; As0; T0; R0];
[t_unc, Y_unc] = ode15s(@(tt,Y) suicide_rhs_uncontrolled(Y, ...
    Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi), ...
    t, Y0);

%% -----------------------------
% Output structure
%% -----------------------------
out.t = t;
out.S = S; out.M = M; out.I = I; out.Af = Af; out.As = As; out.T = Tt; out.R = R;
out.lambda1 = lambda1; out.lambda2 = lambda2; out.lambda3 = lambda3;
out.lambda4 = lambda4; out.lambda5 = lambda5; out.lambda6 = lambda6; out.lambda7 = lambda7;
out.xi1 = xi1; out.xi2 = xi2; out.xi3 = xi3;
out.J = J;
out.stateCost = stateCost;
out.controlCost = controlCost;
out.iter = iter;
out.err = err;
out.strategy = [use_xi1, use_xi2, use_xi3];
out.uncontrolled_t = t_unc;
out.uncontrolled_Y = Y_unc;

%% -----------------------------
% Plots
%% -----------------------------
figure(1);
plot(t, M, 'b', 'LineWidth', 5); hold on;
plot(t_unc, Y_unc(:,2), 'r--', 'LineWidth', 5);
legend('With control','Without control','Location','best');
ylabel('$M$','Interpreter','latex');
xlabel('Time [Years]');
xlim([0 tf]);
box on; grid on; set(gca,'LineWidth',2);

figure(2);
semilogy(t, max(I,1e-12), 'k', 'LineWidth', 5); hold on;
semilogy(t_unc, max(Y_unc(:,3),1e-12), 'g--', 'LineWidth', 5);
legend('With control','Without control','Location','best');
ylabel('$I$','Interpreter','latex');
xlabel('Time [Years]');
xlim([0 tf]);
box on; grid on; set(gca,'LineWidth',2);

figure(3);
semilogy(t, max(As,1e-12), 'c', 'LineWidth', 5); hold on;
semilogy(t_unc, max(Y_unc(:,5),1e-12), 'm--', 'LineWidth', 5);
legend('With control','Without control','Location','best');
ylabel('$A_s$','Interpreter','latex');
xlabel('Time [Years]');
xlim([0 tf]);
box on; grid on; set(gca,'LineWidth',2);

figure(4); hold on;
if use_xi1
    plot(t, xi1, 'b', 'LineWidth', 5);
end
if use_xi2
    plot(t, xi2, 'r', 'LineWidth', 5);
end
if use_xi3
    plot(t, xi3, 'k', 'LineWidth', 5);
end
legendEntries = {};
if use_xi1, legendEntries{end+1} = '$\xi_1$'; end
if use_xi2, legendEntries{end+1} = '$\xi_2$'; end
if use_xi3, legendEntries{end+1} = '$\xi_3$'; end
legend(legendEntries,'Interpreter','latex','Location','best');
ylabel('Control profile');
xlabel('Time [Years]');
ylim([0 1]);
box on; grid on; set(gca,'LineWidth',2);

fprintf('\nState cost   = %.6f\n', stateCost);
fprintf('Control cost = %.6f\n', controlCost);
fprintf('Total cost J = %.6f\n', J);

end

%% ============================================================
% Controlled RHS
%% ============================================================
function dx = suicide_rhs_controlled(x, u, use_xi1, use_xi2, use_xi3, ...
    Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi)

S  = x(1);
M  = x(2);
I  = x(3);
Af = x(4); %#ok<NASGU>
As = x(5);
Tt = x(6);
R  = x(7);

xi1 = u(1); if ~use_xi1, xi1 = 0; end
xi2 = u(2); if ~use_xi2, xi2 = 0; end
xi3 = u(3); if ~use_xi3, xi3 = 0; end

N = S + M + I + As + Tt + R;
if N <= 0
    N = 1e-12;
end

dx = zeros(7,1);
dx(1) = Lambda - (1-xi1)*beta*((I+As)/N)*S - mu*S + phi*omega*R;
dx(2) = (1-xi1)*beta*((I+As)/N)*S - (gamma + (1-xi1)*eta*I + mu)*M + (1-phi)*omega*R;
dx(3) = (gamma + (1-xi1)*eta*I)*M - (alpha + tau*(1+xi2) + mu)*I;
dx(4) = rho*alpha*I;
dx(5) = (1-rho)*alpha*I - (theta*(1+xi3) + mu)*As;
dx(6) = tau*(1+xi2)*I + theta*(1+xi3)*As - (sigma + mu)*Tt;
dx(7) = sigma*Tt - (omega + mu)*R;
end

%% ============================================================
% Uncontrolled RHS
%% ============================================================
function dY = suicide_rhs_uncontrolled(Y, Lambda, mu, beta, eta, gamma, tau, alpha, rho, theta, sigma, omega, phi)

S  = Y(1);
M  = Y(2);
I  = Y(3);
Af = Y(4); %#ok<NASGU>
As = Y(5);
Tt = Y(6);
R  = Y(7);

N = S + M + I + As + Tt + R;
if N <= 0
    N = 1e-12;
end

dY = zeros(7,1);
dY(1) = Lambda - beta*((I+As)/N)*S - mu*S + phi*omega*R;
dY(2) = beta*((I+As)/N)*S - (gamma + eta*I + mu)*M + (1-phi)*omega*R;
dY(3) = (gamma + eta*I)*M - (alpha + tau + mu)*I;
dY(4) = rho*alpha*I;
dY(5) = (1-rho)*alpha*I - (theta + mu)*As;
dY(6) = tau*I + theta*As - (sigma + mu)*Tt;
dY(7) = sigma*Tt - (omega + mu)*R;
end