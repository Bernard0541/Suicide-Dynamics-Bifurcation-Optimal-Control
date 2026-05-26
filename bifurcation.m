%% ============================================================
%  Suicide model bifurcation diagrams (FULL IMPROVED CODE)
%
%  This script generates:
%    1) A forward bifurcation diagram
%    2) A backward bifurcation diagram
%
%  Improvements:
%    - Uses the suicide-model threshold formula for beta(R0)
%    - Computes eta_c from the center-manifold expression
%    - Uses a stronger eta for the backward case
%    - Finds ALL positive endemic equilibria by solving a scalar equation in I*
%    - Classifies branches by the full Jacobian eigenvalues
%    - Uses log-scale for endemic branches only
%    - Marks the critical fold point Rc in the backward diagram
%    - Produces smoother curves with denser grids
%
%  IMPORTANT:
%    The disease-free equilibrium is exactly I*=0, which cannot be shown
%    directly on a log scale. So:
%      - the DFE is plotted on the zero line with plot(...)
%      - endemic branches are plotted with semilogy(...)
%
%% ============================================================

clear; clc; close all;

%% -----------------------------
%  Baseline parameters
%  Replace with calibrated values if available
%% -----------------------------
par.Lambda = 2000;%50;
par.mu     = 0.0092;%0.02;
par.gamma  = 0.300;%0.30;
par.alpha  = 0.7384;%0.20;
par.tau    = 0.3840;%0.15;
par.theta  = 0.5;%0.25;
par.sigma  = 0.400;%0.40;
par.omega  = 0.100;%0.10;
par.rho    = 0.5000;%0.30;
par.phi    = 0.600;%0.60;

%% -----------------------------
%  C-notation
%% -----------------------------
C1 = par.theta + par.mu;
C2 = par.sigma + par.mu;
C3 = par.omega + par.mu;
C4 = par.alpha + par.tau + par.mu;
C5 = 1 - par.rho;
C6 = par.gamma + par.mu;
C7 = 1 - par.phi;

%% -----------------------------
%  beta* from R0 = 1
%% -----------------------------
DenR0 = C4*C1*C6*C2*C3 ...
      - par.gamma*C7*par.omega*par.sigma*(par.tau*C1 + par.theta*C5*par.alpha);

NumR0 = par.gamma*(C1 + C5*par.alpha)*C2*C3;

beta_star = DenR0 / NumR0;

%% -----------------------------
%  eta_c from bifurcation analysis
%% -----------------------------
Psi_tilde = ...
      par.gamma/C4 ...
    + (C5*par.alpha*par.gamma)/(C1*C4) ...
    + (par.gamma^2)/(C4^2) ...
    + (2*C5*par.alpha*par.gamma^2)/(C1*C4^2) ...
    + (C5^2*par.alpha^2*par.gamma^2)/(C1^2*C4^2) ...
    + (par.gamma^2*(par.tau*C1 + par.theta*C5*par.alpha))/(C1*C2*C4^2) ...
    + (C5*par.alpha*par.gamma^2*(par.tau*C1 + par.theta*C5*par.alpha))/(C1^2*C2*C4^2) ...
    + (par.sigma*par.gamma^2*(par.tau*C1 + par.theta*C5*par.alpha))/(C1*C2*C3*C4^2) ...
    + (C5*par.alpha*par.sigma*par.gamma^2*(par.tau*C1 + par.theta*C5*par.alpha))/(C1^2*C2*C3*C4^2);

eta_c = beta_star * C4 * Psi_tilde / par.Lambda;

fprintf('beta* = %.10f\n', beta_star);
fprintf('eta_c = %.10f\n', eta_c);

%% -----------------------------
%  Choose eta values
%% -----------------------------
eta_forward  = 0.50 * eta_c;
eta_backward = 100 * eta_c;   % strong backward choice for visible subcritical branch

fprintf('eta_forward  = %.10f\n', eta_forward);
fprintf('eta_backward = %.10f\n', eta_backward);

%% -----------------------------
%  R0 grids
%% -----------------------------
R0_values_global   = linspace(0, 2.5, 2500);
R0_values_backward = linspace(0.88, 2.5, 2500);

%% -----------------------------
%  Compute forward branches
%% -----------------------------
fprintf('Computing forward branches...\n');
dataF = compute_branches_scalar(R0_values_global, par, eta_forward);

%% -----------------------------
%  Compute backward branches
%% -----------------------------
fprintf('Computing backward branches...\n');
dataB = compute_branches_scalar(R0_values_backward, par, eta_backward);

%% -----------------------------
%  Compute Rc for backward case
%% -----------------------------
Rc = NaN;
Ic = NaN;
if ~isempty(dataB.R0_unstable)
    [Rc, idxRc] = max(dataB.R0_unstable);
    Ic = dataB.I_unstable(idxRc);
end
fprintf('Estimated Rc = %.10f\n', Rc);

%% -----------------------------
%  Visualization settings
%% -----------------------------
I_floor = 1e-4;   % positive floor for log visualization only

%% -----------------------------
%  Main figure
%% -----------------------------
figure('Color','w','Position',[80 80 1450 560]);

%% ============================================================
%  PANEL 1: Forward bifurcation
%% ============================================================
subplot(1,2,1); 
hold on; box on;

% DFE
plot(R0_values_global(R0_values_global<=1), zeros(sum(R0_values_global<=1),1), ...
    'k-', 'LineWidth', 2.0);
plot(R0_values_global(R0_values_global>1), zeros(sum(R0_values_global>1),1), ...
    'k--', 'LineWidth', 2.0);

% Sort and plot endemic branches
[RFs, IFs] = sort_branch(dataF.R0_stable, dataF.I_stable);
[RFu, IFu] = sort_branch(dataF.R0_unstable, dataF.I_unstable);

if ~isempty(RFs)
    semilogy(RFs, max(IFs, I_floor), 'b-', 'LineWidth', 4);
end
if ~isempty(RFu)
    semilogy(RFu, max(IFu, I_floor), 'r--', 'LineWidth', 4);
end

xline(1,'k:','LineWidth',1.4);

xlabel('$\mathcal{R}_0$','Interpreter','latex');
ylabel('$I^*$','Interpreter','latex');
title(sprintf('Forward bifurcation ($\\eta = %.4g < \\eta_c$)', eta_forward), ...
    'Interpreter','latex');

% legend({'DFE stable','DFE unstable','SPE stable','SPE unstable'}, ...
%     'Location','northwest');

set(gca,'YScale','log');
xlim([0.8 2.5]);
%ylim([I_floor 300]);
ylim([15 4000]);
grid on;
%hold off

%text(0.06, 1.5*I_floor, '$I^*=0$', 'Interpreter','latex', 'FontSize', 11);

%% ============================================================
%  PANEL 2: Backward bifurcation
%% ============================================================
%figure('Color','w','Position',[80 80 1450 560]);
subplot(1,2,2); 

hold on; box on;

% DFE
plot(R0_values_backward(R0_values_backward<=1), zeros(sum(R0_values_backward<=1),1), ...
    'k-', 'LineWidth', 2.0);
plot(R0_values_backward(R0_values_backward>1), zeros(sum(R0_values_backward>1),1), ...
    'k--', 'LineWidth', 2.0);

% Sort and plot endemic branches
[RBs, IBs] = sort_branch(dataB.R0_stable, dataB.I_stable);
[RBu, IBu] = sort_branch(dataB.R0_unstable, dataB.I_unstable);

if ~isempty(RBs)
    semilogy(RBs, max(IBs, I_floor), 'b-', 'LineWidth', 4);
end
if ~isempty(RBu)
    semilogy(RBu, max(IBu, I_floor), 'r--', 'LineWidth', 4);
end

xline(1,'k:','LineWidth',2);
%xline(0.933157,'k:','LineWidth',2);

% Mark Rc
% if ~isnan(Rc)
%     plot(Rc, max(Ic, I_floor), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
%     xline(Rc, 'm--', 'LineWidth', 1.4);
%     text(Rc + 0.002, max(Ic, I_floor)*1.25, ...
%         sprintf('$\\mathcal{R}_c = %.4f$', Rc), ...
%         'Interpreter','latex', 'Color','m', 'FontSize', 11);
% end

xlabel('$\mathcal{R}_0$','Interpreter','latex');
ylabel('$I^*$','Interpreter','latex');
title(sprintf('Backward bifurcation ($\\eta = %.4g > \\eta_c$)', eta_backward), ...
    'Interpreter','latex');

% if ~isnan(Rc)
%     legend({'DFE stable','DFE unstable','SPE stable','SPE unstable','$\mathcal{R}_c$'}, ...
%         'Interpreter','latex', 'Location','northwest');
% else
%     legend({'DFE stable','DFE unstable','SPE stable','SPE unstable'}, ...
%         'Interpreter','latex', 'Location','northwest');
% end

set(gca,'YScale','log');
xlim([0.88 1.05]);
%xlim([0 2]);
%ylim([I_floor 80]);
ylim([0.1 2000]);
grid on;
%hold off

%text(0.885, 1.5*I_floor, '$I^*=0$', 'Interpreter','latex', 'FontSize', 11);

%% -----------------------------
%  Optional separate zoom figure
%% -----------------------------
% figure('Color','w','Position',[180 120 700 560]); hold on; box on;
% 
% plot(R0_values_backward(R0_values_backward<=1), zeros(sum(R0_values_backward<=1),1), ...
%     'k-', 'LineWidth', 2.0);
% plot(R0_values_backward(R0_values_backward>1), zeros(sum(R0_values_backward>1),1), ...
%     'k--', 'LineWidth', 2.0);
% 
% if ~isempty(RBs)
%     semilogy(RBs, max(IBs, I_floor), 'b-', 'LineWidth', 2.6);
% end
% if ~isempty(RBu)
%     semilogy(RBu, max(IBu, I_floor), 'r--', 'LineWidth', 2.4);
% end
% 
% xline(1,'k:','LineWidth',1.4);
% if ~isnan(Rc)
%     plot(Rc, max(Ic, I_floor), 'ko', 'MarkerFaceColor','y', 'MarkerSize',8);
%     xline(Rc, 'm--', 'LineWidth', 1.5);
%     text(Rc + 0.0015, max(Ic, I_floor)*1.30, ...
%         sprintf('$\\mathcal{R}_c = %.5f$', Rc), ...
%         'Interpreter','latex', 'Color','m', 'FontSize', 12);
% end
% 
% xlabel('$\mathcal{R}_0$','Interpreter','latex');
% ylabel('$I^*$','Interpreter','latex');
% title('Backward bifurcation (zoomed view)','Interpreter','latex');
% 
% set(gca,'YScale','log');
% xlim([0.90 1.01]);
% %ylim([I_floor 50]);
% ylim([0.001 50]);
% grid on;
% 
% if ~isnan(Rc)
%     legend({'DFE stable','DFE unstable','SPE stable','SPE unstable','$\mathcal{R}_c$'}, ...
%         'Interpreter','latex', 'Location','northwest');
% else
%     legend({'DFE stable','DFE unstable','SPE stable','SPE unstable'}, ...
%         'Interpreter','latex', 'Location','northwest');
% end

%text(0.902, 1.5*I_floor, '$I^*=0$', 'Interpreter','latex', 'FontSize', 11);

%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================

function data = compute_branches_scalar(R0_values, par, eta_val)

    stable_R0   = [];
    stable_I    = [];
    unstable_R0 = [];
    unstable_I  = [];

    for kk = 1:length(R0_values)

        R0 = R0_values(kk);
        rootsI = find_positive_I_roots_given_R0(par, eta_val, R0);

        for j = 1:length(rootsI)
            Ieq = rootsI(j);
            xeq = equilibrium_from_I(Ieq, par, eta_val, R0);

            if isempty(xeq), continue; end
            if any(~isfinite(xeq)), continue; end
            if any(xeq < -1e-8), continue; end

            J = suicide_jacobian(xeq, par, eta_val, R0);
            ev = eig(J);

            if max(real(ev)) < -1e-7
                stable_R0(end+1,1) = R0;       %#ok<AGROW>
                stable_I(end+1,1)  = xeq(3);   %#ok<AGROW>
            else
                unstable_R0(end+1,1) = R0;     %#ok<AGROW>
                unstable_I(end+1,1)  = xeq(3); %#ok<AGROW>
            end
        end
    end

    data.R0_stable   = stable_R0;
    data.I_stable    = stable_I;
    data.R0_unstable = unstable_R0;
    data.I_unstable  = unstable_I;
end

function rootsI = find_positive_I_roots_given_R0(par, eta_val, R0)
    rootsI = [];

    if par.rho * par.alpha > 0
        Imax = 0.9999 * par.Lambda / (par.rho * par.alpha);
    else
        Imax = 500;
    end

    Imin  = 1e-12;
    Igrid = linspace(Imin, Imax, 50000);
    Hvals = arrayfun(@(I) H_scalar(I, par, eta_val, R0), Igrid);

    mask = isfinite(Hvals);
    Igrid = Igrid(mask);
    Hvals = Hvals(mask);

    % sign-change roots
    for k = 1:length(Igrid)-1
        a = Igrid(k);
        b = Igrid(k+1);
        fa = Hvals(k);
        fb = Hvals(k+1);

        root = NaN;
        if fa == 0
            root = a;
        elseif fa * fb < 0
            try
                root = fzero(@(I) H_scalar(I, par, eta_val, R0), [a b]);
            catch
                root = NaN;
            end
        end

        if isfinite(root) && root > 1e-8
            if isempty(rootsI) || all(abs(rootsI - root) > 1e-7)
                rootsI(end+1,1) = root; %#ok<AGROW>
            end
        end
    end

    % near-tangent roots
    AbsH = abs(Hvals);
    for k = 2:length(Igrid)-1
        if AbsH(k) < AbsH(k-1) && AbsH(k) < AbsH(k+1) && AbsH(k) < 1e-7
            a = max(Igrid(max(1,k-30)), Imin);
            b = min(Igrid(min(length(Igrid),k+30)), Imax);
            try
                root = fzero(@(I) H_scalar(I, par, eta_val, R0), [a b]);
            catch
                root = NaN;
            end

            if isfinite(root) && root > 1e-8
                if isempty(rootsI) || all(abs(rootsI - root) > 1e-7)
                    rootsI(end+1,1) = root; %#ok<AGROW>
                end
            end
        end
    end

    rootsI = sort(rootsI);
end

function H = H_scalar(I, par, eta_val, R0)

    C1 = par.theta + par.mu;
    C2 = par.sigma + par.mu;
    C3 = par.omega + par.mu;
    C4 = par.alpha + par.tau + par.mu;
    C5 = 1 - par.rho;
    C6 = par.gamma + par.mu;
    C7 = 1 - par.phi;

    DenR0 = C4*C1*C6*C2*C3 ...
          - par.gamma*C7*par.omega*par.sigma*(par.tau*C1 + par.theta*C5*par.alpha);

    NumR0 = par.gamma*(C1 + C5*par.alpha)*C2*C3;
    beta  = R0 * DenR0 / NumR0;

    if I <= 0
        H = NaN;
        return;
    end

    a = C5 * par.alpha / C1;            % As = a I
    t = (par.tau + par.theta*a) / C2;   % T = t I
    r = par.sigma * t / C3;             % R = r I
    q = 1 + a + t + r;

    % From N' = Lambda - mu*N - rho*alpha*I = 0
    N = (par.Lambda - par.rho * par.alpha * I) / par.mu;
    if N <= 0
        H = NaN;
        return;
    end

    % From I' = 0
    denomM = par.gamma + eta_val * I;
    if denomM <= 0
        H = NaN;
        return;
    end
    M = C4 * I / denomM;

    % From S' = 0
    denomS = par.Lambda + I * (beta * (1 + a) - par.rho * par.alpha);
    if abs(denomS) < 1e-14
        H = NaN;
        return;
    end
    S = N * (par.Lambda + par.phi * par.omega * r * I) / denomS;

    % Closure
    H = N - S - M - q * I;
end

function xeq = equilibrium_from_I(I, par, eta_val, R0)

    C1 = par.theta + par.mu;
    C2 = par.sigma + par.mu;
    C3 = par.omega + par.mu;
    C4 = par.alpha + par.tau + par.mu;
    C5 = 1 - par.rho;
    C6 = par.gamma + par.mu;
    C7 = 1 - par.phi;

    DenR0 = C4*C1*C6*C2*C3 ...
          - par.gamma*C7*par.omega*par.sigma*(par.tau*C1 + par.theta*C5*par.alpha);

    NumR0 = par.gamma*(C1 + C5*par.alpha)*C2*C3;
    beta  = R0 * DenR0 / NumR0;

    a = C5 * par.alpha / C1;
    t = (par.tau + par.theta * a) / C2;
    r = par.sigma * t / C3;

    N = (par.Lambda - par.rho * par.alpha * I) / par.mu;
    if N <= 0
        xeq = [];
        return;
    end

    M  = C4 * I / (par.gamma + eta_val * I);
    As = a * I;
    T  = t * I;
    R  = r * I;
    S  = N - M - I - As - T - R;

    xeq = [S; M; I; As; T; R];
end

function J = suicide_jacobian(x, par, eta_val, R0)

    S  = x(1);
    M  = x(2);
    I  = x(3);
    As = x(4);
    T  = x(5);
    R  = x(6);

    C1 = par.theta + par.mu;
    C2 = par.sigma + par.mu;
    C3 = par.omega + par.mu;
    C4 = par.alpha + par.tau + par.mu;
    C5 = 1 - par.rho;
    C6 = par.gamma + par.mu;
    C7 = 1 - par.phi;

    DenR0 = C4*C1*C6*C2*C3 ...
          - par.gamma*C7*par.omega*par.sigma*(par.tau*C1 + par.theta*C5*par.alpha);

    NumR0 = par.gamma*(C1 + C5*par.alpha)*C2*C3;
    beta  = R0 * DenR0 / NumR0;

    N = S + M + I + As + T + R;
    if N <= 0
        N = 1e-12;
    end

    P = I + As;

    % h = beta*S*(I+As)/N
    dh_dS  = beta * P * (N - S) / N^2;
    dh_dM  = -beta * S * P / N^2;
    dh_dI  = beta * S * (N - P) / N^2;
    dh_dAs = beta * S * (N - P) / N^2;
    dh_dT  = -beta * S * P / N^2;
    dh_dR  = -beta * S * P / N^2;

    J = zeros(6,6);

    % f1
    J(1,1) = -dh_dS - par.mu;
    J(1,2) = -dh_dM;
    J(1,3) = -dh_dI;
    J(1,4) = -dh_dAs;
    J(1,5) = -dh_dT;
    J(1,6) = -dh_dR + par.phi * par.omega;

    % f2
    J(2,1) = dh_dS;
    J(2,2) = dh_dM - (par.gamma + eta_val * I + par.mu);
    J(2,3) = dh_dI - eta_val * M;
    J(2,4) = dh_dAs;
    J(2,5) = dh_dT;
    J(2,6) = dh_dR + (1 - par.phi) * par.omega;

    % f3
    J(3,1) = 0;
    J(3,2) = par.gamma + eta_val * I;
    J(3,3) = eta_val * M - (par.alpha + par.tau + par.mu);
    J(3,4) = 0;
    J(3,5) = 0;
    J(3,6) = 0;

    % f4
    J(4,1) = 0;
    J(4,2) = 0;
    J(4,3) = (1 - par.rho) * par.alpha;
    J(4,4) = -(par.theta + par.mu);
    J(4,5) = 0;
    J(4,6) = 0;

    % f5
    J(5,1) = 0;
    J(5,2) = 0;
    J(5,3) = par.tau;
    J(5,4) = par.theta;
    J(5,5) = -(par.sigma + par.mu);
    J(5,6) = 0;

    % f6
    J(6,1) = 0;
    J(6,2) = 0;
    J(6,3) = 0;
    J(6,4) = 0;
    J(6,5) = par.sigma;
    J(6,6) = -(par.omega + par.mu);
end

function [Rsorted, Isorted] = sort_branch(R, I)
    if isempty(R)
        Rsorted = [];
        Isorted = [];
        return;
    end

    A = [R(:), I(:)];
    A = unique(round(A, 10), 'rows');
    A = sortrows(A, [1 2]);

    Rsorted = A(:,1);
    Isorted = A(:,2);
end