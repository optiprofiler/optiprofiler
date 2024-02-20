function runAllInOne(solvers, labels, problem_names, varargin)
    %RUNBENCHMARK plots performance profiles and data profiles.

    % Preprocess the solvers.
    % TODO
    % Preprocess the labels.
    % TODO

    % Set default values for optional arguments.
    myCluster = parcluster('local');
    nb_cores = myCluster.NumWorkers;
    profile_options = struct(ProfileOptionKey.N_JOBS.value, nb_cores, ProfileOptionKey.BENCHMARK_ID.value, '.');
    problem_options = struct();

    if mod(length(varargin), 2) ~= 0
        error('Options should be provided as name-value pairs.');
    end
    for i = 1:2:length(varargin)
        key = varargin{i};
        value = varargin{i + 1};

        validProblemOptionKeys = {enumeration('ProblemOptionKey').value};
        validProfileOptionKeys = {enumeration('ProfileOptionKey').value};

        if ismember(key, validProblemOptionKeys)
            problem_options.(key) = value;
        elseif ismember(key, validProfileOptionKeys)
            profile_options.(key) = value;
        else
            error("Unknown option: %s", key);
        end
    end

    % Judge whether profile_options.n_jobs is a integer between 1 and nb_cores.
    if isfield(profile_options, ProfileOptionKey.N_JOBS.value)
        if ~isnumeric(profile_options.(ProfileOptionKey.N_JOBS.value))
            error("profile_options.n_jobs should be a integer.");
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) < 1
            profile_options.(ProfileOptionKey.N_JOBS.value) = 1;
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) > nb_cores
            profile_options.(ProfileOptionKey.N_JOBS.value) = nb_cores;
        else
            profile_options.(ProfileOptionKey.N_JOBS.value) = round(profile_options.(ProfileOptionKey.N_JOBS.value));
        end
    end

    % Judge whether profile_options.benchmark_id is a string.
    if ~isstring(profile_options.(ProfileOptionKey.BENCHMARK_ID.value)) && ~ischar(profile_options.(ProfileOptionKey.BENCHMARK_ID.value))
        error("profile_options.benchmark_id should be a string.");
    end

    % Build features.
    features = cell(6, 1);
    features{1} = Feature('plain');
    features{2} = Feature('noisy');
    features{3} = Feature('randomize_x0');
    features{4} = Feature('regularized');
    features{5} = Feature('tough');
    features{6} = Feature('truncated');

    % Paths to the results.
    timestamp = datestr(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH-mm-SSZ'), 'yyyy-mm-ddTHH-MM-SSZ');
    full_path = mfilename('fullpath');
    [folder_path, ~, ~] = fileparts(full_path);
    root_path = fileparts(folder_path);
    path_out = fullfile(root_path, 'out', 'one-for-all', profile_options.(ProfileOptionKey.BENCHMARK_ID.value), timestamp);
    if ~exist(path_out, 'dir')
        mkdir(path_out);
    end
    pdf_summary_hist = fullfile(path_out, 'summary_hist.pdf');
    pdf_summary_ret = fullfile(path_out, 'summary_ret.pdf');
    
    % Store the names of the problems.
    path_txt = fullfile(path_out, 'problems.txt');
    fid = fopen(path_txt, 'w');
    if fid == -1
        error("Cannot open the file %s.", path_txt);
    end
    fprintf(fid, '%s\n', sort(string(problem_names)));
    fclose(fid);

    max_eval_factor = 500;
    % Run the benchmark for each feature.
    for i_feature = 1:length(features)
        feature = features{i_feature};
        fprintf("INFO: Starting the computation of the %s profiles.\n", feature.name);
        
        % Solve all the problems.
        [fun_histories, maxcv_histories, fun_ret, maxcv_ret, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions] = solveAllProblems(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);
        merit_histories = computeMeritValues(fun_histories, maxcv_histories);
        merit_ret = computeMeritValues(fun_ret, maxcv_ret);
        merit_init = computeMeritValues(fun_init, maxcv_init);

        % Determine the least merit value for each problem.
        merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
        if ismember(feature.name, {FeatureName.NOISY.value, FeatureName.TOUGH.value, FeatureName.TRUNCATED.value})
            feature_plain = Feature(FeatureName.PLAIN.value);
            fprintf("INFO: Starting the computation of the plain profiles.\n");
            [fun_histories_plain, maxcv_histories_plain] = solveAllProblems(problem_names, problem_options, solvers, labels, feature_plain, max_eval_factor, profile_options);
            merit_histories_plain = computeMeritValues(fun_histories_plain, maxcv_histories_plain);
            merit_min_plain = min(min(min(merit_histories_plain, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
            merit_min = min(merit_min, merit_min_plain, 'omitnan');
        end

        % Set the default values for plotting.
        fprintf("Creating results.\n");
        set(groot, 'DefaultLineLineWidth', 1);
        set(groot, 'DefaultAxesFontSize', 12);
        set(groot, 'DefaultAxesFontName', 'Arial');
        set(groot, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
        set(groot, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});
        
        [n_problems, n_solvers, n_runs, max_eval] = size(merit_histories);

        % Create the figure for the summary.
        defaultFigurePosition = get(0, 'DefaultFigurePosition');
        default_width = defaultFigurePosition(3);
        default_height = defaultFigurePosition(4);
        fig_summary_hist = figure('Position', [defaultFigurePosition(1:2), 10 * default_width, 2 * default_height], 'visible', 'off');
        fig_summary_ret = figure('Position', [defaultFigurePosition(1:2), 10 * default_width, 2 * default_height], 'visible', 'off');
        t_summary_hist = tiledlayout(fig_summary_hist, 2, 10, 'Padding', 'compact', 'TileSpacing', 'compact');
        t_summary_ret = tiledlayout(fig_summary_ret, 2, 10, 'Padding', 'compact', 'TileSpacing', 'compact');
        for i = 1:20
            axs_summary_hist(i) = nexttile(t_summary_hist);
            axs_summary_ret(i) = nexttile(t_summary_ret);
        end

        tolerances = 10.^(-1:-1:-10);
        for i_profile = 1:10
            tolerance = tolerances(i_profile);
            tolerance_label = ['$\tau = 10^{', int2str(log10(tolerance)), '}$'];

            work_hist = NaN(n_problems, n_solvers, n_runs);
            work_ret = NaN(n_problems, n_solvers, n_runs);
            for i_problem = 1:n_problems
                for i_solver = 1:n_solvers
                    for i_run = 1:n_runs
                        if isfinite(merit_min(i_problem))
                            threshold = max(tolerance * merit_init(i_problem) + (1 - tolerance) * merit_min(i_problem), merit_min(i_problem));
                        else
                            threshold = -Inf;
                        end
                        if min(merit_histories(i_problem, i_solver, i_run, :), [], 'omitnan') <= threshold
                            work_hist(i_problem, i_solver, i_run) = find(merit_histories(i_problem, i_solver, i_run, :) <= threshold, 1, 'first');
                        end
                        if merit_ret(i_problem, i_solver, i_run) <= threshold
                            work_ret(i_problem, i_solver, i_run) = n_eval(i_problem, i_solver, i_run);
                        end
                    end
                end
            end

            % Draw the performance and data profiles.
            [x_perf_hist, y_perf_hist, ratio_max_perf_hist, x_data_hist, y_data_hist, ratio_max_data_hist] = getExtendedPerformancesDataProfileAxes(work_hist, problem_dimensions);
            [x_perf_ret, y_perf_ret, ratio_max_perf_ret, x_data_ret, y_data_ret, ratio_max_data_ret] = getExtendedPerformancesDataProfileAxes(work_ret, problem_dimensions);
            drawPerformanceDataProfiles(axs_summary_hist(i_profile), x_perf_hist, y_perf_hist, labels);
            drawPerformanceDataProfiles(axs_summary_ret(i_profile), x_perf_ret, y_perf_ret, labels);
            drawPerformanceDataProfiles(axs_summary_hist(10 + i_profile), x_data_hist, y_data_hist, labels);
            drawPerformanceDataProfiles(axs_summary_ret(10 + i_profile), x_data_ret, y_data_ret, labels);

            % Set x-axis limits.
            set(axs_summary_hist(i_profile), 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
            set(axs_summary_ret(i_profile), 'XLim', [0.0, 1.1 * ratio_max_perf_ret]);
            set(axs_summary_hist(10 + i_profile), 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
            set(axs_summary_ret(10 + i_profile), 'XLim', [0.0, 1.1 * ratio_max_data_ret]);

            % Modify x-axis ticks labels of the performance profiles.
            xtl_summary_1 = arrayfun(@(x) ['$', '2^{', num2str(x), '}$'], get(axs_summary_hist(i_profile), 'XTick'), 'UniformOutput', false);
            set(axs_summary_hist(i_profile), 'XTickLabel', xtl_summary_1, 'TickLabelInterpreter', 'latex');
            xtl_summary_2 = arrayfun(@(x) ['$', '2^{', num2str(x), '}$'], get(axs_summary_ret(i_profile), 'XTick'), 'UniformOutput', false);
            set(axs_summary_ret(i_profile), 'XTickLabel', xtl_summary_2, 'TickLabelInterpreter', 'latex');

            % Set x-axis labels.
            xlabel(axs_summary_hist(i_profile), 'Performance ratio', 'Interpreter', 'latex');
            xlabel(axs_summary_ret(i_profile), 'Performance ratio', 'Interpreter', 'latex');
            xlabel(axs_summary_hist(10 + i_profile), 'Number of simplex gradients', 'Interpreter', 'latex');
            xlabel(axs_summary_ret(10 + i_profile), 'Number of simplex gradients', 'Interpreter', 'latex');

            % Set y-axis labels.
            ylabel(axs_summary_hist(i_profile), ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            ylabel(axs_summary_ret(i_profile), ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            ylabel(axs_summary_hist(10 + i_profile), ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            ylabel(axs_summary_ret(10 + i_profile), ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        end

        if i_feature == 1
            exportgraphics(fig_summary_hist, pdf_summary_hist, 'ContentType', 'vector');
            exportgraphics(fig_summary_ret, pdf_summary_ret, 'ContentType', 'vector');
        else
            exportgraphics(fig_summary_hist, pdf_summary_hist, 'ContentType', 'vector', 'Append', true);
            exportgraphics(fig_summary_ret, pdf_summary_ret, 'ContentType', 'vector', 'Append', true);
        end

        close(fig_summary_hist);
        close(fig_summary_ret);

    end

end


function merit_values = computeMeritValues(fun_values, maxcv_values)
    is_nearly_feasible = maxcv_values <= 1e-12;
    is_very_infeasible = maxcv_values >= 1e-6;
    is_undecided = ~is_nearly_feasible & ~is_very_infeasible;
    merit_values = NaN(size(fun_values));
    merit_values(is_nearly_feasible) = fun_values(is_nearly_feasible);
    merit_values(is_very_infeasible) = Inf;
    merit_values(is_undecided) = fun_values(is_undecided) + 1e8 * maxcv_values(is_undecided);
end