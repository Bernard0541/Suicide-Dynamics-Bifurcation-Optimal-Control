function cea = suicide_cost_effectiveness_analysis
%% ============================================================
% Cost-effectiveness analysis for the suicide optimal control model
%
% Strategies:
%   A = (xi1, xi2, xi3)
%   B = (xi1, xi2)
%   C = (xi1, xi3)
%   D = (xi2, xi3)
%
% Effectiveness:
%   Cases averted = uncontrolled burden - controlled burden
%   burden = integral of (M + I + As) over [0, tf]
%
% Cost:
%   Control implementation cost only
%   C(xi) = int_0^tf 0.5*(C1*xi1^2 + C2*xi2^2 + C3*xi3^2) dt
%
% Figures:
%   Figure 1: Total cost and ACER only for A, B, C, D,
%             arranged in increasing order of total cost
%   Figure 2 onward: Total cost, ACER, and ICER
%                    arranged in increasing order of total cost and
%                    eliminating the least cost-effective strategy iteratively.
%
% NOTE:
%   To avoid NaN or Inf in the ICER plots, ACER is used as a replacement:
%   (i) for the first strategy in the total-cost ordering, and
%   (ii) whenever incremental effectiveness is non-positive.
%% ============================================================

clearvars; clc; close all;

%% -----------------------------
% Parameters
%% -----------------------------
par.tf     = 49.0;
par.nSteps = 5000;
par.t      = linspace(0, par.tf, par.nSteps+1);
par.h      = par.tf / par.nSteps;

par.Lambda = 2000;
par.mu     = 0.0092;
par.beta   = 0.9000;
par.eta    = 1.0000;
par.gamma  = 0.3000;
par.alpha  = 0.7384;
par.rho    = 0.5000;
par.theta  = 0.5000;
par.sigma  = 0.4000;
par.omega  = 0.1000;
par.phi    = 0.6000;
par.tau    = 0.3840;   % psi in LaTeX = tau in MATLAB code

%% -----------------------------
% Cost weights in objective
%% -----------------------------
par.B1 = 120;
par.B2 = 150;
par.B3 = 100;

par.C1 = 100;
par.C2 = 100;
par.C3 = 100;

%% -----------------------------
% Initial conditions
%% -----------------------------
par.S0  = 99980;
par.M0  = 1;
par.I0  = 4;
par.Af0 = 0;
par.As0 = 2;
par.T0  = 1;
par.R0  = 1;

%% -----------------------------
% Numerical settings
%% -----------------------------
par.tol        = 1e-6;
par.maxIter    = 500;
par.relaxParam = 0.1;

%% -----------------------------
% Uncontrolled baseline
%% -----------------------------
base = run_uncontrolled(par);
fprintf('\nUncontrolled burden = %.6f\n', base.burden);

%% -----------------------------
% Strategies
%% -----------------------------
strategies = {
    'A', true,  true,  true;
    'B', true,  true,  false;
    'C', true,  false, true;
    'D', false, true,  true
    };

nStrat = size(strategies,1);

names        = strings(nStrat,1);
costs        = zeros(nStrat,1);
burdens      = zeros(nStrat,1);
casesAverted = zeros(nStrat,1);
ACER         = nan(nStrat,1);

allResults = cell(nStrat,1);

%% -----------------------------
% Run all strategies
%% -----------------------------
for s = 1:nStrat
    name    = strategies{s,1};
    use_xi1 = strategies{s,2};
    use_xi2 = strategies{s,3};
    use_xi3 = strategies{s,4};

    fprintf('\nRunning Strategy %s ...\n', name);

    out = run_controlled_strategy(par, use_xi1, use_xi2, use_xi3);

    names(s)        = name;
    costs(s)        = out.controlCost;
    burdens(s)      = out.burden;
    casesAverted(s) = base.burden - out.burden;

    if casesAverted(s) > 0
        ACER(s) = costs(s) / casesAverted(s);
    else
        ACER(s) = 0;
    end

    allResults{s} = out;

    fprintf('  Burden         = %.6f\n', out.burden);
    fprintf('  Cases averted  = %.6f\n', casesAverted(s));
    fprintf('  Control cost   = %.6f\n', costs(s));
    fprintf('  ACER           = %.6f\n', ACER(s));
end

%% -----------------------------
% Tables
% Arrange in increasing order of total cost
%% -----------------------------
[~, idxTable] = sort(costs, 'ascend');

ACER_Table = table(names(idxTable), burdens(idxTable), casesAverted(idxTable), ...
    costs(idxTable), ACER(idxTable), ...
    'VariableNames', {'Strategy','TotalBurden','CasesAverted','ControlCost','ACER'});

disp(' ');
disp('================ ACER TABLE ================');
disp(ACER_Table);

%% -----------------------------
% Figure 1: Total cost and ACER only
% Arrange strategies in increasing order of total cost
%% -----------------------------
namesFig1 = names(idxTable);
costsFig1 = costs(idxTable);
ACERFig1  = ACER(idxTable);

cats1 = categorical(namesFig1);
cats1 = reordercats(cats1, cellstr(namesFig1));

figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 900 400]);

subplot(1,2,1);
bar(cats1, costsFig1, 'FaceColor', [0.2 0.5 0.8]);
ylabel('Total cost');
title('Total cost for Strategies A-D');
box on;

subplot(1,2,2);
bar(cats1, ACERFig1, 'FaceColor', [0.3 0.7 0.4]);
ylabel('ACER');
title('ACER for Strategies A-D');
box on;

%% -----------------------------
% Best strategy by ACER
%% -----------------------------
[bestACER, idxBestACER] = min(ACER(ACER>0));
idxCandidates = find(ACER == bestACER, 1, 'first');
bestACERStrategy = names(idxCandidates);

fprintf('\nMost cost-effective by ACER: %s (ACER = %.6f)\n', bestACERStrategy, bestACER);

%% -----------------------------
% ICER iterative elimination
% Before each ICER computation, sort by increasing total cost
%% -----------------------------
currentNames = names;
currentCosts = costs;
currentAvert = casesAverted;

roundCount = 2;
eliminationLog = strings(0,1);
ICER_Rounds = {};

while numel(currentNames) > 1
    [curSorted, curICER, curACER] = compute_current_icer(currentNames, currentCosts, currentAvert);

    ICER_Rounds{end+1} = table(curSorted.names, curSorted.costs, curSorted.casesAverted, ...
        curACER, curICER, ...
        'VariableNames', {'Strategy','ControlCost','CasesAverted','ACER','ICER'}); %#ok<AGROW>

    make_icer_figure(roundCount, curSorted.names, curSorted.costs, curSorted.casesAverted, curACER, curICER);

    idxRemove = choose_strategy_to_remove(curSorted.names, curSorted.costs, curSorted.casesAverted, curICER);

    strategyRemoved = curSorted.names(idxRemove);
    eliminationLog(end+1,1) = strategyRemoved; %#ok<AGROW>

    fprintf('\nFigure %d analysis:\n', roundCount);
    fprintf('  Strategies arranged in increasing order of total cost:\n');
    for i = 1:numel(curSorted.names)
        fprintf('    %s : cost = %.6f, cases averted = %.6f, ACER = %.6f, ICER = %.6f\n', ...
            curSorted.names(i), curSorted.costs(i), curSorted.casesAverted(i), curACER(i), curICER(i));
    end
    fprintf('  Removing least cost-effective strategy: %s\n', strategyRemoved);

    keepMask = curSorted.names ~= strategyRemoved;
    currentNames = curSorted.names(keepMask);
    currentCosts = curSorted.costs(keepMask);
    currentAvert = curSorted.casesAverted(keepMask);

    roundCount = roundCount + 1;
end

finalICERStrategy = currentNames(1);
fprintf('\nFinal strategy selected by iterative ICER elimination: %s\n', finalICERStrategy);

%% -----------------------------
% Final summary
%% -----------------------------
Final_Summary = table(names, costs, casesAverted, ACER, ...
    'VariableNames', {'Strategy','ControlCost','CasesAverted','ACER'});

disp(' ');
disp('================ FINAL SUMMARY ================');
disp(Final_Summary);

%% -----------------------------
% Return results
%% -----------------------------
cea.base = base;
cea.results = allResults;
cea.ACER_Table = ACER_Table;
cea.bestACERStrategy = bestACERStrategy;
cea.bestACER = bestACER;
cea.finalICERStrategy = finalICERStrategy;
cea.eliminationLog = eliminationLog;
cea.ICER_Rounds = ICER_Rounds;

end

%% ============================================================
% Compute sequential ICER on current set
% Sort by increasing total cost before ICER calculation
% Replace NaN/Inf with ACER values
%% ============================================================
function [sortedData, ICER, ACER_local] = compute_current_icer(names, costs, casesAverted)

    [sortedCosts, idx] = sort(costs, 'ascend');
    namesSorted  = names(idx);
    avertSorted  = casesAverted(idx);

    n = numel(namesSorted);
    ICER = zeros(n,1);
    ACER_local = zeros(n,1);

    % ACER for all current strategies
    for i = 1:n
        if avertSorted(i) > 0
            ACER_local(i) = sortedCosts(i) / avertSorted(i);
        else
            ACER_local(i) = 0;
        end
    end

    % First strategy: replace ICER with ACER
    ICER(1) = ACER_local(1);

    % Remaining strategies
    for i = 2:n
        dCost = sortedCosts(i) - sortedCosts(i-1);
        dEff  = avertSorted(i) - avertSorted(i-1);

        if dEff > 0
            ICER(i) = dCost / dEff;
        else
            % Replace problematic ICER with ACER
            ICER(i) = ACER_local(i);
        end
    end

    sortedData.names = namesSorted;
    sortedData.costs = sortedCosts;
    sortedData.casesAverted = avertSorted;

end

%% ============================================================
% Choose strategy to remove in current ICER round
%% ============================================================
function idxRemove = choose_strategy_to_remove(names, costs, casesAverted, ICER)

n = numel(names);

% Strong dominance:
% remove a strategy if another strategy is cheaper and more effective
for i = 1:n
    for j = 1:n
        if i ~= j
            if (costs(i) >= costs(j)) && (casesAverted(i) <= casesAverted(j)) && ...
               ((costs(i) > costs(j)) || (casesAverted(i) < casesAverted(j)))
                idxRemove = i;
                return;
            end
        end
    end
end

% Otherwise remove the strategy with the largest ICER
tempICER = ICER;
tempICER(1) = -inf;  % cheapest comparator

[~, idxRemove] = max(tempICER);

end

%% ============================================================
% Make ICER figure for current round
%% ============================================================
function make_icer_figure(figNum, names, costs, casesAverted, ACER_local, ICER)

figure(figNum); clf;
set(gcf, 'Color', 'w', 'Position', [120 120 1350 400]);

cats = categorical(names);
cats = reordercats(cats, cellstr(names));

subplot(1,3,1);
bar(cats, costs, 'FaceColor', [0.2 0.5 0.8]);
ylabel('Total cost');
title('Total cost');
box on;

subplot(1,3,2);
bar(cats, ACER_local, 'FaceColor', [0.3 0.7 0.4]);
ylabel('ACER');
title('ACER');
box on;

subplot(1,3,3);
bar(cats, ICER, 'FaceColor', [0.8 0.4 0.4]);
ylabel('ICER');
title('ICER');
box on;

end

%% ============================================================
% Controlled strategy runner
%% ============================================================
function out = run_controlled_strategy(par, use_xi1, use_xi2, use_xi3)

t      = par.t;
h      = par.h;
nSteps = par.nSteps;

S  = zeros(1,nSteps+1);
M  = zeros(1,nSteps+1);
I  = zeros(1,nSteps+1);
Af = zeros(1,nSteps+1);
As = zeros(1,nSteps+1);
Tt = zeros(1,nSteps+1);
R  = zeros(1,nSteps+1);

S(1)=par.S0; M(1)=par.M0; I(1)=par.I0; Af(1)=par.Af0;
As(1)=par.As0; Tt(1)=par.T0; R(1)=par.R0;

lambda1 = zeros(1,nSteps+1);
lambda2 = zeros(1,nSteps+1);
lambda3 = zeros(1,nSteps+1);
lambda4 = zeros(1,nSteps+1);
lambda5 = zeros(1,nSteps+1);
lambda6 = zeros(1,nSteps+1);
lambda7 = zeros(1,nSteps+1);

xi1 = 0.5*ones(1,nSteps+1);
xi2 = 0.5*ones(1,nSteps+1);
xi3 = 0.5*ones(1,nSteps+1);

if ~use_xi1, xi1(:)=0; end
if ~use_xi2, xi2(:)=0; end
if ~use_xi3, xi3(:)=0; end

iter = 0;
err  = inf;

while err > par.tol && iter < par.maxIter
    iter = iter + 1;

    oldS = S; oldM = M; oldI = I; oldAf = Af; oldAs = As; oldTt = Tt; oldR = R;
    oldlambda1 = lambda1; oldlambda2 = lambda2; oldlambda3 = lambda3;
    oldlambda4 = lambda4; oldlambda5 = lambda5; oldlambda6 = lambda6; oldlambda7 = lambda7;
    oldxi1 = xi1; oldxi2 = xi2; oldxi3 = xi3;

    for k = 1:nSteps
        xk = [S(k); M(k); I(k); Af(k); As(k); Tt(k); R(k)];
        uk = [xi1(k); xi2(k); xi3(k)];

        k1 = suicide_rhs_controlled(xk, uk, use_xi1, use_xi2, use_xi3, par);
        k2 = suicide_rhs_controlled(xk + 0.5*h*k1, uk, use_xi1, use_xi2, use_xi3, par);
        k3 = suicide_rhs_controlled(xk + 0.5*h*k2, uk, use_xi1, use_xi2, use_xi3, par);
        k4 = suicide_rhs_controlled(xk + h*k3, uk, use_xi1, use_xi2, use_xi3, par);

        xnext = xk + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        xnext = max(xnext,0);

        S(k+1)=xnext(1); M(k+1)=xnext(2); I(k+1)=xnext(3); Af(k+1)=xnext(4);
        As(k+1)=xnext(5); Tt(k+1)=xnext(6); R(k+1)=xnext(7);
    end

    if any(~isfinite([S M I Af As Tt R]))
        error('State variables became non-finite.');
    end

    lambda1(end)=0; lambda2(end)=0; lambda3(end)=0; lambda4(end)=0;
    lambda5(end)=0; lambda6(end)=0; lambda7(end)=0;

    for k = nSteps:-1:1
        Nk = S(k+1)+M(k+1)+I(k+1)+As(k+1)+Tt(k+1)+R(k+1);
        if Nk <= 0
            Nk = 1e-12;
        end

        xi1k = xi1(k+1); if ~use_xi1, xi1k = 0; end
        xi2k = xi2(k+1); if ~use_xi2, xi2k = 0; end
        xi3k = xi3(k+1); if ~use_xi3, xi3k = 0; end

        dl1 = par.beta*(lambda1(k+1)-lambda2(k+1))*(1-xi1k)*(I(k+1)+As(k+1))/Nk ...
              + par.mu*lambda1(k+1);

        dl2 = (lambda2(k+1)-lambda3(k+1))*(par.gamma + (1-xi1k)*par.eta*I(k+1)) ...
              + par.mu*lambda2(k+1) - par.B1;

        dl3 = par.alpha*(lambda3(k+1)-par.rho*lambda4(k+1)) ...
              - par.alpha*(1-par.rho)*lambda5(k+1) ...
              + (lambda3(k+1)-lambda6(k+1))*(1+xi2k)*par.tau ...
              + par.mu*lambda3(k+1) - par.B2 ...
              + (lambda2(k+1)-lambda3(k+1))*(1-xi1k)*par.eta*M(k+1) ...
              + par.beta*(lambda1(k+1)-lambda2(k+1))*(1-xi1k)*(I(k+1)+As(k+1))/Nk;

        dl4 = 0;

        dl5 = par.beta*(lambda1(k+1)-lambda2(k+1))*(1-xi1k)*S(k+1)/Nk ...
              + (lambda5(k+1)-lambda6(k+1))*(1+xi3k)*par.theta ...
              + par.mu*lambda5(k+1) - par.B3;

        dl6 = par.sigma*(lambda6(k+1)-lambda7(k+1)) + par.mu*lambda6(k+1);

        dl7 = (par.omega+par.mu)*lambda7(k+1) ...
              - (1-par.phi)*par.omega*lambda2(k+1) - par.phi*par.omega*lambda1(k+1);

        lambda1(k)=lambda1(k+1)-h*dl1;
        lambda2(k)=lambda2(k+1)-h*dl2;
        lambda3(k)=lambda3(k+1)-h*dl3;
        lambda4(k)=lambda4(k+1)-h*dl4;
        lambda5(k)=lambda5(k+1)-h*dl5;
        lambda6(k)=lambda6(k+1)-h*dl6;
        lambda7(k)=lambda7(k+1)-h*dl7;
    end

    if any(~isfinite([lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7]))
        error('Adjoint variables became non-finite.');
    end

    Nall = S + M + I + As + Tt + R;
    Nall(Nall <= 0) = 1e-12;

    xi1_proj = (par.eta*(lambda3-lambda2).*I.*M.*Nall + par.beta*(lambda2-lambda1).*(I+As).*S) ./ (par.C1*Nall);
    xi2_proj = (par.tau*(lambda3-lambda6).*I) ./ par.C2;
    xi3_proj = (par.theta*(lambda5-lambda6).*As) ./ par.C3;

    xi1_proj = min(1,max(0,xi1_proj));
    xi2_proj = min(1,max(0,xi2_proj));
    xi3_proj = min(1,max(0,xi3_proj));

    if ~use_xi1, xi1_proj(:)=0; end
    if ~use_xi2, xi2_proj(:)=0; end
    if ~use_xi3, xi3_proj(:)=0; end

    xi1 = par.relaxParam*xi1_proj + (1-par.relaxParam)*oldxi1;
    xi2 = par.relaxParam*xi2_proj + (1-par.relaxParam)*oldxi2;
    xi3 = par.relaxParam*xi3_proj + (1-par.relaxParam)*oldxi3;

    if ~use_xi1, xi1(:)=0; end
    if ~use_xi2, xi2(:)=0; end
    if ~use_xi3, xi3(:)=0; end

    err = max([
        norm(xi1-oldxi1, inf), norm(xi2-oldxi2, inf), norm(xi3-oldxi3, inf), ...
        norm(S-oldS, inf), norm(M-oldM, inf), norm(I-oldI, inf), ...
        norm(Af-oldAf, inf), norm(As-oldAs, inf), norm(Tt-oldTt, inf), norm(R-oldR, inf), ...
        norm(lambda1-oldlambda1, inf), norm(lambda2-oldlambda2, inf), norm(lambda3-oldlambda3, inf), ...
        norm(lambda4-oldlambda4, inf), norm(lambda5-oldlambda5, inf), norm(lambda6-oldlambda6, inf), ...
        norm(lambda7-oldlambda7, inf)
        ]);
end

if iter == par.maxIter && err > par.tol
    warning('Maximum iterations reached before convergence.');
end

burden = trapz(t, M + I + As);
controlCost = trapz(t, 0.5*par.C1*xi1.^2 + 0.5*par.C2*xi2.^2 + 0.5*par.C3*xi3.^2);

out.t = t;
out.S = S; out.M = M; out.I = I; out.Af = Af; out.As = As; out.T = Tt; out.R = R;
out.xi1 = xi1; out.xi2 = xi2; out.xi3 = xi3;
out.controlCost = controlCost;
out.burden = burden;
out.iter = iter;
out.err = err;

end

%% ============================================================
% Uncontrolled baseline
%% ============================================================
function base = run_uncontrolled(par)

Y0 = [par.S0; par.M0; par.I0; par.Af0; par.As0; par.T0; par.R0];

[t_unc, Y_unc] = ode15s(@(tt,Y) suicide_rhs_uncontrolled(Y, par), par.t, Y0);

burden_unc = trapz(t_unc, Y_unc(:,2) + Y_unc(:,3) + Y_unc(:,5));

base.t = t_unc;
base.Y = Y_unc;
base.burden = burden_unc;

end

%% ============================================================
% Controlled RHS
%% ============================================================
function dx = suicide_rhs_controlled(x, u, use_xi1, use_xi2, use_xi3, par)

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
dx(1) = par.Lambda - (1-xi1)*par.beta*((I+As)/N)*S - par.mu*S + par.phi*par.omega*R;
dx(2) = (1-xi1)*par.beta*((I+As)/N)*S - (par.gamma + (1-xi1)*par.eta*I + par.mu)*M + (1-par.phi)*par.omega*R;
dx(3) = (par.gamma + (1-xi1)*par.eta*I)*M - (par.alpha + par.tau*(1+xi2) + par.mu)*I;
dx(4) = par.rho*par.alpha*I;
dx(5) = (1-par.rho)*par.alpha*I - (par.theta*(1+xi3) + par.mu)*As;
dx(6) = par.tau*(1+xi2)*I + par.theta*(1+xi3)*As - (par.sigma + par.mu)*Tt;
dx(7) = par.sigma*Tt - (par.omega + par.mu)*R;

end

%% ============================================================
% Uncontrolled RHS
%% ============================================================
function dY = suicide_rhs_uncontrolled(Y, par)

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
dY(1) = par.Lambda - par.beta*((I+As)/N)*S - par.mu*S + par.phi*par.omega*R;
dY(2) = par.beta*((I+As)/N)*S - (par.gamma + par.eta*I + par.mu)*M + (1-par.phi)*par.omega*R;
dY(3) = (par.gamma + par.eta*I)*M - (par.alpha + par.tau + par.mu)*I;
dY(4) = par.rho*par.alpha*I;
dY(5) = (1-par.rho)*par.alpha*I - (par.theta + par.mu)*As;
dY(6) = par.tau*I + par.theta*As - (par.sigma + par.mu)*Tt;
dY(7) = par.sigma*Tt - (par.omega + par.mu)*R;

end