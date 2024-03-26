function drawLogRatioProfiles(ax, work, labels)
%DRAWLOGRATIOPROFILES draws log-ratio profiles.

    [n_problems, n_solvers, n_runs] = size(work);
    work_flat = reshape(permute(work, [1, 3, 2]), n_problems * n_runs, n_solvers);
    log_ratio = NaN(n_problems * n_runs, 1);
    log_ratio_finite = isfinite(work_flat(:, 1)) & isfinite(work_flat(:, 2));
    log_ratio(log_ratio_finite) = log2(work_flat(log_ratio_finite, 1)) - log2(work_flat(log_ratio_finite, 2));
    ratio_max = max(max(abs(log_ratio(log_ratio_finite)), [], 'all'), eps);
    if isempty(ratio_max)
        ratio_max = eps;
    end
    log_ratio(isnan(work_flat(:, 1)) & isfinite(work_flat(:, 2))) = 1.1 * ratio_max;
    log_ratio(isfinite(work_flat(:, 1)) & isnan(work_flat(:, 2))) = -1.1 * ratio_max;
    log_ratio(isnan(work_flat(:, 1)) & isnan(work_flat(:, 2))) = 0.0;
    log_ratio = sort(log_ratio);

    x = (1:(n_problems * n_runs))';
    bar(ax, x(log_ratio < 0), log_ratio(log_ratio < 0), 'LineStyle', 'None', 'LineWidth', 1.5);
    hold(ax, 'on');
    bar(ax, x(log_ratio > 0), log_ratio(log_ratio > 0), 'LineStyle', 'None', 'LineWidth', 1.5);
    text(ax, (n_problems * n_runs + 1) / 2, -ratio_max, labels{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
    text(ax, (n_problems * n_runs + 1) / 2, ratio_max, labels{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);
    xticks(ax, []);
    xlim(ax, [0.5, n_problems * n_runs + 0.5]);
    yyaxis(ax, 'right');
    set(ax, "YTickLabel", []);
    yyaxis(ax, 'left');
    ylim(ax, [-1.1 * ratio_max, 1.1 * ratio_max]);
    linkprop([ax.XAxis; ax.YAxis],'color');
    linkprop([ax.YAxis(1), ax.YAxis(2)],{'Limits','TickValues'});
    box(ax,'on');
    hold(ax, 'off');

end