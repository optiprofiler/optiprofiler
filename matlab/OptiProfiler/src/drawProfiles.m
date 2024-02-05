function [fig_summary, fig_perf_hist, fig_perf_out, fig_data_hist, fig_data_out, fig_log_ratio_hist, fig_log_ratio_out] = drawProfiles(work_hist, work_out, problem_dimensions, labels, tolerance_label)

    n_solvers = size(work_hist, 2);

    % Create the figure for the summary.
    defaultFigurePosition = get(0, 'DefaultFigurePosition');
    default_width = defaultFigurePosition(3);
    default_height = defaultFigurePosition(4);
    if n_solvers > 2
        fig_summary = figure('Position', [defaultFigurePosition(1:2), 2 * default_width, 2 * default_height]);
        for i = 1:4
            axs_summary(i) = subplot(2, 2, i);
        end
    else
        fig_summary = figure('Position', [defaultFigurePosition(1:2), 2 * default_width, 3 * default_height]);
        for i = 1:6
            axs_summary(i) = subplot(3, 2, i);
        end
    end

    % Create the individual figures.
    fig_perf_hist = figure;
    ax_perf_hist = axes('Parent', fig_perf_hist);
    fig_perf_out = figure;
    ax_perf_out = axes('Parent', fig_perf_out);
    fig_data_hist = figure;
    ax_data_hist = axes('Parent', fig_data_hist);
    fig_data_out = figure;
    ax_data_out = axes('Parent', fig_data_out);
    if n_solvers > 2
        fig_log_ratio_hist = [];
        ax_log_ratio_hist = [];
        fig_log_ratio_out = [];
        ax_log_ratio_out = [];
    else
        fig_log_ratio_hist = figure;
        ax_log_ratio_hist = axes('Parent', fig_log_ratio_hist);

        fig_log_ratio_out = figure;
        ax_log_ratio_out = axes('Parent', fig_log_ratio_out);
    end
    

end