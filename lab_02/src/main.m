X = [-8.47, -7.45, -7.12, -8.30, -8.15, -6.01, -5.20, -7.38, -6.76, -9.18, -6.00, -8.08, -7.96, -8.34, -6.82, -8.46, -8.07, -7.04, -7.24, -8.16, -8.20, -8.27, -7.79, -7.37, -7.02, -7.13, -6.99, -7.65, -8.18, -6.71, -8.41, -6.71, -7.04, -9.15, -7.74, -10.11, -8.20, -7.07, -7.63, -8.99, -6.62, -6.23, -7.13, -6.41, -7.06, -7.72, -8.44, -8.85, -8.02, -6.98, -6.08, -7.20, -7.48, -7.82, -9.19, -8.31, -7.95, -7.97, -6.66, -6.59, -9.10, -7.87, -9.02, -8.77, -7.62, -9.44, -8.05, -7.60, -7.33, -6.94, -8.51, -7.39, -6.44, -8.88, -8.21, -7.66, -6.91, -8.39, -7.37, -7.26, -6.04, -7.58, -7.28, -7.02, -7.10, -7.33, -8.63, -8.21, -7.12, -8.11, -9.03, -8.11, -8.79, -9.22, -7.32, -5.97, -7.26, -6.39, -7.64, -8.38, -7.67, -7.70, -7.70, -8.95, -6.25, -8.09, -7.85, -8.10, -7.73, -6.78, -7.78, -8.20, -8.88, -8.51, -7.45, -7.14, -6.63, -7.38, -7.72, -6.25];
gamma = 0.9;

% 1-2
[muhat, muci] = my_normfit_mu(X, 1 - gamma);
[s2hat, s2ci] = my_normfit_s2(X, 1 - gamma);
fprintf("mu = %f\nmu_bottom = %f\nmu_top = %f\n", muhat, muci(1), muci(2));
fprintf("s2 = %f\ns2_bottom = %f\ns2_top = %f\n", s2hat, s2ci(1), s2ci(2));
% 3
process_mu(X, gamma, muhat);
process_s2(X, gamma, s2hat);


function [muhat, muci] = normfit_mu(X, alpha)
    [muhat, ~, muci, ~] = normfit(X, alpha);
end

function [s2hat, s2ci] = normfit_s2(X, alpha)
    [~, sigmahat, ~, sigmaci] = normfit(X, alpha);
    s2hat = sigmahat ^ 2;
    s2ci = sigmaci .^ 2;
end

function [muhat, muci] = my_normfit_mu(X, alpha)
    muhat = mean(X);
    s = std(X);
    gamma = 1 - alpha;
    n = length(X);
    mu_bottom = muhat + s * tinv((1 - gamma) / 2, n - 1) / sqrt(n);
    mu_top = muhat + s * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
    muci = [mu_bottom, mu_top];
end

function [s2hat, s2ci] = my_normfit_s2(X, alpha)
    s2hat = var(X);
    gamma = 1 - alpha;
    n = length(X);
    s2_top = (n - 1) * s2hat / chi2inv((1 - gamma) / 2, n - 1);
    s2_bottom = (n - 1) * s2hat / chi2inv((1 + gamma) / 2, n - 1);
    s2ci = [s2_bottom, s2_top];
end

function process_parameter(start, X, gamma, est, fit, line_legend, est_legend, top_legend, bottom_legend)
    N = length(X);
    figure;
    hold on;
    grid on;
    plot([1, N], [est, est]);
    ests = [];
    cis_bottom = [];
    cis_top = [];
    for n = start:N
        [est, cis] = fit(X(1:n), 1 - gamma);
        ests = [ests, est];
        cis_bottom = [cis_bottom, cis(1)];
        cis_top = [cis_top, cis(2)];
    end
    plot(start:N, ests);
    plot(start:N, cis_bottom);
    plot(start:N, cis_top);
    l = legend(line_legend, est_legend, top_legend, bottom_legend);
    set(l, 'Interpreter', 'latex', 'fontsize', 18);
    hold off;
end

function process_mu(X, gamma, muhat)
    process_parameter(1, X, gamma, muhat, @my_normfit_mu, '$\hat\mu(\vec x_N)$', '$\hat\mu(\vec x_n)$', '$\underline\mu(\vec x_n)$', '$\overline\mu(\vec x_n)$');
end

function process_s2(X, gamma, S2)
    process_parameter(10, X, gamma, S2, @my_normfit_s2, '$\hat\sigma^2(\vec x_N)$', '$\hat\sigma^2(\vec x_n)$', '$\underline\sigma^2(\vec x_n)$', '$\overline\sigma^2(\vec x_n)$');
end
