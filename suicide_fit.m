function results = suicide_fit(dataFile)
% Improved parameter estimation for the suicide model
%
% Fits annual suicide death rate per 100,000 using yearly increments
% of the cumulative death state A_f.
%
% Additions:
%   1) Weighted least squares
%   2) Residual bootstrap confidence band (95% CI)
%   3) CI plotted on the Observed vs Fitted figure
%
% Model states:
%   x = [S, M, I, A_s, T, R, A_f]
%
% Data file should contain columns:
%   Year, Rate, Count
%
% Example:
%   results = suicide_fit('suicide_death.csv');

clearvars; clc

    if nargin < 1
        dataFile = 'suicide_death.csv';
    end

    %% -----------------------------
    %  Read data
    %% -----------------------------
    T = readtable(dataFile);

    years = T.Year(:);
    yobs  = T.Rate(:);   % annual suicide death rate per 100,000

    nObs = numel(yobs);
    tEdges = 0:nObs;     % 25 points for 24 yearly increments

    %% -----------------------------
    %  Initial conditions (fixed)
    %  Interpreted on a per-100,000 basis
    %% -----------------------------
    x0 = [99980; 12; 4; 2; 1; 1; 0];  % [S, M, I, As, T, R, Af]

    %% -----------------------------
    %  Fixed parameters
    %% -----------------------------
    fixed.Lambda = 2000;
    fixed.mu     = 0.020;
    fixed.gamma  = 0.300;
    fixed.sigma  = 0.400;
    fixed.omega  = 0.100;
    fixed.phi    = 0.600;

    %% -----------------------------
    %  Parameters to estimate
    %% -----------------------------
    parNames = {'beta','eta','alpha','tau','rho','theta'};
    p0 = [0.25, 0.001, 0.20, 0.15, 0.30, 0.25];

    lb = [0, 0,    0, 0,    0, 0];
    ub = [1.0,  1.0,  1.0,  5.0,  1.0, 1.0];

    %% -----------------------------
    %  Weighted least squares weights
    %  Using inverse-scale weights to reduce heteroscedasticity
    %% -----------------------------
    w = 1 ./ max(yobs, 0.5);   % column vector

    %% -----------------------------
    %  Multi-start least squares
    %% -----------------------------
    nStarts = 20;
    bestP   = [];
    bestRes = [];
    bestSSE = inf;
    bestResnorm = inf;

    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxFunctionEvaluations', 2e4, ...
        'MaxIterations', 2e3, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);

    rng(1); % reproducibility

    for k = 1:nStarts
        if k == 1
            pInit = p0;
        else
            pInit = lb + rand(size(lb)).*(ub-lb);
        end

        try
            [pHat, resnorm, residual] = lsqnonlin( ...
                @(p) residuals_fun(p, parNames, fixed, x0, tEdges, yobs, w), ...
                pInit, lb, ub, opts);

            % Compute unweighted SSE for reporting
            parsTmp = fixed;
            for i = 1:numel(parNames)
                parsTmp.(parNames{i}) = pHat(i);
            end
            [~, Xtmp] = simulate_model(parsTmp, x0, tEdges);
            yhatTmp = diff(Xtmp(:,7));
            SSE = sum((yhatTmp - yobs).^2);

            if resnorm < bestResnorm
                bestResnorm = resnorm;
                bestSSE     = SSE;
                bestP       = pHat;
                bestRes     = residual;
            end
        catch ME
            warning('Start %d failed: %s', k, ME.message);
        end
    end

    if isempty(bestP)
        error('All optimization starts failed.');
    end

    %% -----------------------------
    %  Final fitted parameter structure
    %% -----------------------------
    pars = fixed;
    for i = 1:numel(parNames)
        pars.(parNames{i}) = bestP(i);
    end

    %% -----------------------------
    %  Final simulation
    %% -----------------------------
    [~, X] = simulate_model(pars, x0, tEdges);

    % Annual predicted deaths per 100,000 from cumulative Af increments
    yhat = diff(X(:,7));

    %% -----------------------------
    %  Summary statistics
    %% -----------------------------
    SSE  = sum((yhat - yobs).^2);
    RMSE = sqrt(mean((yhat - yobs).^2));
    MAE  = mean(abs(yhat - yobs));

    kPar = numel(bestP);
    n    = numel(yobs);
    AIC  = n*log(SSE/n) + 2*kPar;
    BIC  = n*log(SSE/n) + kPar*log(n);

    %% -----------------------------
    %  Bootstrap 95% confidence band
    %  Residual bootstrap on the data scale
    %% -----------------------------
    nBoot = 1000; % raise to 500+ for final publication runs
    yboot = nan(nObs, nBoot);
    pboot = nan(kPar, nBoot);

    fprintf('\nRunning bootstrap for 95%% CI with %d replicates...\n', nBoot);

    rawRes = yobs - yhat;
    centeredRes = rawRes - mean(rawRes);

    bootOpts = optimoptions('lsqnonlin', ...
        'Display','off', ...
        'MaxFunctionEvaluations', 1e4, ...
        'MaxIterations', 1e3, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);

    for b = 1:nBoot
        resampled = centeredRes(randi(nObs, nObs, 1));
        ybootObs = yhat + resampled;
        ybootObs = max(ybootObs, 0);  % keep rates nonnegative

        try
            pBoot = lsqnonlin( ...
                @(p) residuals_fun(p, parNames, fixed, x0, tEdges, ybootObs, w), ...
                bestP, lb, ub, bootOpts);

            parsBoot = fixed;
            for i = 1:numel(parNames)
                parsBoot.(parNames{i}) = pBoot(i);
            end

            [~, XBoot] = simulate_model(parsBoot, x0, tEdges);
            yboot(:,b) = diff(XBoot(:,7));
            pboot(:,b) = pBoot(:);
        catch
            % leave failed replicate as NaN
        end
    end

    % Remove failed bootstrap fits
    good = all(isfinite(yboot),1);
    yboot = yboot(:,good);
    pboot = pboot(:,good);

    if isempty(yboot)
        warning('All bootstrap runs failed. CI not available.');
        CI_low  = nan(size(yhat));
        CI_high = nan(size(yhat));
        ystd    = nan(size(yhat));
    else
        CI_low  = prctile(yboot',  2.5)';
        CI_high = prctile(yboot', 97.5)';
        ystd    = std(yboot, 0, 2);
    end

    %% -----------------------------
    %  Display results
    %% -----------------------------
    fprintf('\nEstimated parameters:\n');
    for i = 1:numel(parNames)
        fprintf('  %-6s = %.6f\n', parNames{i}, bestP(i));
    end
    fprintf('  SSE   = %.6f\n', SSE);
    fprintf('  RMSE  = %.6f\n', RMSE);
    fprintf('  MAE   = %.6f\n', MAE);
    fprintf('  AIC   = %.6f\n', AIC);
    fprintf('  BIC   = %.6f\n', BIC);
    fprintf('  Bootstrap successful runs = %d / %d\n\n', size(yboot,2), nBoot);

    %% -----------------------------
    %  Plot observed vs fitted with 95% CI
    %% -----------------------------
    figure('Color','w'); hold on;

    if all(isfinite(CI_low)) && all(isfinite(CI_high))
        fill([years; flipud(years)], [CI_low; flipud(CI_high)], ...
             [1.0 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.7);
    end

    plot(years, yhat, 'r-', 'LineWidth',5);
    plot(years, yobs, 'ko', 'MarkerFaceColor','b', 'MarkerSize',7);

    xlabel('Year');
    ylabel('Suicide death rate per 100,000');
    %title('Observed vs fitted annual suicide death rate with 95% CI');
    xlim([2000 2023]); xticks(2000:2:2023);

    if all(isfinite(CI_low)) && all(isfinite(CI_high))
        legend('95% CI','Cumulative death ($A_f$)','Data','Location','best');
    else
        legend('Fitted','Observed','Location','best');
    end
    box on;

    %% -----------------------------
    %  Plot state trajectories
    %% -----------------------------
    % figure('Color','w');
    % plot(tEdges, X(:,1), 'b-', 'LineWidth',2); hold on;
    % plot(tEdges, X(:,2), 'r-', 'LineWidth',2);
    % plot(tEdges, X(:,3), 'c-', 'LineWidth',2);
    % plot(tEdges, X(:,4), 'm-', 'LineWidth',2);
    % plot(tEdges, X(:,5), 'g-', 'LineWidth',2);
    % plot(tEdges, X(:,6), 'y-', 'LineWidth',2);
    % plot(tEdges, X(:,7), 'k-', 'LineWidth',2);
    % xlabel('Time (years since 2000)');
    % ylabel('Population per 100,000');
    % legend('S','M','I','A_s','T','R','A_f','Location','best');
    % title('Model trajectories under fitted parameters');
    % box on;

    %% -----------------------------
    %  Return results
    %% -----------------------------
    results = struct();
    results.parNames = parNames;
    results.pHat     = bestP;
    results.fixed    = fixed;
    results.SSE      = SSE;
    results.RMSE     = RMSE;
    results.MAE      = MAE;
    results.AIC      = AIC;
    results.BIC      = BIC;
    results.weights  = w;
    results.years    = years;
    results.yobs     = yobs;
    results.yhat     = yhat;
    results.CI_low   = CI_low;
    results.CI_high  = CI_high;
    results.ystd     = ystd;
    results.yboot    = yboot;
    results.pboot    = pboot;
    results.t        = tEdges;
    results.X        = X;
    results.residual = bestRes;
end

%% ============================================================
% Weighted residual function
%% ============================================================
function r = residuals_fun(p, parNames, fixed, x0, tEdges, yobs, w)
    pars = fixed;
    for i = 1:numel(parNames)
        pars.(parNames{i}) = p(i);
    end

    [~, X] = simulate_model(pars, x0, tEdges);

    % Annual deaths per 100,000 from increments of cumulative Af
    yhat = diff(X(:,7));

    if any(~isfinite(yhat)) || any(yhat < 0)
        r = 1e6 * ones(size(yobs));
        return;
    end

    % Weighted least squares residuals
    r = (yhat - yobs) .* w;
end

%% ============================================================
% ODE simulation
%% ============================================================
function [tSol, XSol] = simulate_model(pars, x0, tEval)

    opts = odeset('NonNegative',1:7,'RelTol',1e-8,'AbsTol',1e-10);

    [tSol, XSol] = ode15s(@(t,x) suicide_ode(t,x,pars), tEval, x0, opts);

    if size(XSol,1) ~= numel(tEval)
        XSol = interp1(tSol, XSol, tEval, 'pchip');
        tSol = tEval(:);
    end
end

%% ============================================================
% Suicide ODE system
%% ============================================================
function dx = suicide_ode(~, x, p)
    % States: x = [S, M, I, As, T, R, Af]

    S  = x(1);
    M  = x(2);
    I  = x(3);
    As = x(4);
    T  = x(5);
    R  = x(6);

    N = S + M + I + As + T + R;
    if N <= 0
        N = 1e-12;
    end

    lambda = (I + As) / N;

    dx = zeros(7,1);
    dx(1) = p.Lambda - p.beta*lambda*S - p.mu*S + p.phi*p.omega*R;
    dx(2) = p.beta*lambda*S - (p.gamma + p.eta*I + p.mu)*M + (1-p.phi)*p.omega*R;
    dx(3) = (p.gamma + p.eta*I)*M - (p.alpha + p.tau + p.mu)*I;
    dx(4) = (1-p.rho)*p.alpha*I - (p.theta + p.mu)*As;
    dx(5) = p.tau*I + p.theta*As - (p.sigma + p.mu)*T;
    dx(6) = p.sigma*T - (p.omega + p.mu)*R;
    dx(7) = p.rho*p.alpha*I;  % cumulative suicide deaths
end