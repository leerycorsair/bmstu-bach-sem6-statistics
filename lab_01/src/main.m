function main()
    function myhist(X, bins, counts, R, m)
        n = length(X);
        delta = R / m;
        middles = zeros(1, m);
        xx = zeros(1, m);

        for i = 1:m
            xx(i) = counts(i) / (n * delta);
        end

        for i = 1:m
            middles(i) = bins(i) + (delta/2);
        end

        fprintf("    delta: %f\n", delta);
        fprintf("    hist coords:\n");

        for i = 1:m
            fprintf("    [%d] : %f %f\n", i,middles(i), xx(i));
        end

        fprintf("    hist_s = %f\n", sum(xx) * delta);
        set(gca, "xlim", [min(bins) - 1, max(bins) + 1]);
        bar(middles, xx,1, "facecolor", "b", "edgecolor", "w");

        X_range = m_min-delta:(sigma / 100):m_max+delta;
        X_pdf = normpdf(X_range, mu, sigma);
        plot(X_range, X_pdf, "r");
    end

    function mycdf(X, bins, counts)
        n = length(X);
        xx = zeros(1, m + 3);

        bins = [(min(bins) - 0.5) bins (max(bins) + 1)];
        counts = [0 counts 0];

        m = m + 2;
        acc = 0;

        for i = 2:m
            acc = acc + counts(i);
            xx(i) = acc / n;
        end

        xx(m + 1) = 1;

        X_n = (min(X) - 0.5):(sigma / 100):(max(X) + 1.5);
        X_cdf = normcdf(X_n, mu, sigma);
        plot(X_n, X_cdf, "r");

        for i = 2:m
            fprintf("x = %f : F(x) = %f\n", bins(i), xx(i));
        end

        set(gca, "ylim", [0, 1.1]);
        stairs(gca, bins, xx, "b");
    end

    X = [2, 2, 2, 2, 2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,-8.47, -7.45, -7.12, -8.30, -8.15, -6.01, -5.20, -7.38, -6.76, -9.18, -6.00, -8.08, -7.96, -8.34, -6.82, -8.46, -8.07, -7.04, -7.24, -8.16, -8.20, -8.27, -7.79, -7.37, -7.02, -7.13, -6.99, -7.65, -8.18, -6.71, -8.41, -6.71, -7.04, -9.15, -7.74, -10.11, -8.20, -7.07, -7.63, -8.99, -6.62, -6.23, -7.13, -6.41, -7.06, -7.72, -8.44, -8.85, -8.02, -6.98, -6.08, -7.20, -7.48, -7.82, -9.19, -8.31, -7.95, -7.97, -6.66, -6.59, -9.10, -7.87, -9.02, -8.77, -7.62, -9.44, -8.05, -7.60, -7.33, -6.94, -8.51, -7.39, -6.44, -8.88, -8.21, -7.66, -6.91, -8.39, -7.37, -7.26, -6.04, -7.58, -7.28, -7.02, -7.10, -7.33, -8.63, -8.21, -7.12, -8.11, -9.03, -8.11, -8.79, -9.22, -7.32, -5.97, -7.26, -6.39, -7.64, -8.38, -7.67, -7.70, -7.70, -8.95, -6.25, -8.09, -7.85, -8.10, -7.73, -6.78, -7.78, -8.20, -8.88, -8.51, -7.45, -7.14, -6.63, -7.38, -7.72, -6.25];
    m_max = max(X);
    m_min = min(X);
    fprintf("(a) \n");
    fprintf("    max_value = %f\n", m_max);
    fprintf("    min_value = %f\n", m_min);
    fprintf("----------------------------------------\n");

    r = m_max - m_min;
    fprintf("(b) \n");
    fprintf("    r = %f\n", r);
    fprintf("----------------------------------------\n");

    n = length(X);
    mu = sum(X) / n;
    s_2 = sum((X - mu).^2) / (n - 1);
    sigma = sqrt(s_2);

    fprintf("(c) \n");
    fprintf("    mu = %f\n", mu);
    fprintf("    s_2 = %f\n", s_2);
    fprintf("----------------------------------------\n");

    m = floor(log2(n)) + 2;
    bins = [];
    cur = m_min;

    for i = 1:(m + 1)
        bins(i) = cur;
        cur = cur + r / m;
    end

    eps = 1e-6;
    counts = [];
    j = 1;

    for i = 1:(m - 1)
        cur_count = 0;
        for j = 1:n
            if (bins(i) < X(j) || abs(bins(i) - X(j)) < eps) && X(j) < bins(i + 1)
                cur_count = cur_count + 1;
            end
        end
        counts(i) = cur_count;
    end

    cur_count = 0;

    for j = 1:n
        if (bins(m) < X(j) || abs(bins(m) - X(j)) < eps) && (X(j) < bins(m + 1) || abs(bins(m + 1) - X(j)) < eps)
            cur_count = cur_count + 1;
        end
    end

    counts(m) = cur_count;

    fprintf("(d)\n");

    for i = 1:(m)
        fprintf("    [%f : %f) - %d values.\n", bins(i), bins(i + 1), counts(i));
    end

    fprintf("----------------------------------------\n");

    fprintf("(e) \n");

    figure;
    hold on;
    grid on;
    myhist(X, bins, counts, r, m);
    xlabel('X')
    ylabel('P')
    hold off;

    fprintf("----------------------------------------\n");

    fprintf("(f) \n");

    figure;
    hold on;
    grid on;
    mycdf(X, bins, counts);
    xlabel('X')
    ylabel('F')
    hold off;
end