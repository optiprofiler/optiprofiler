function runBenchmark(solvers, labels, problem_names, feature_name, varargin)
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
    feature_options = struct();

    if mod(length(varargin), 2) ~= 0
        error('Options should be provided as name-value pairs.');
    end
    for i = 1:2:length(varargin)
        key = varargin{i};
        value = varargin{i + 1};

        validFeatureOptionKeys = {enumeration('FeatureOptionKey').value};
        validProblemOptionKeys = {enumeration('ProblemOptionKey').value};
        validProfileOptionKeys = {enumeration('ProfileOptionKey').value};

        if ismember(key, validFeatureOptionKeys)
            feature_options.(key) = value;
        elseif ismember(key, validProblemOptionKeys)
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

    % Build feature.
    feature = Feature(feature_name, feature_options);
    fprintf("INFO: Starting the computation of the %s profiles.\n", feature.name);

    % Solve all the problems.
    max_eval_factor = 500;
    [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions] = solveAll(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);
    merit_histories = computeMeritValues(fun_histories, maxcv_histories);
    merit_out = computeMeritValues(fun_out, maxcv_out);
    merit_init = computeMeritValues(fun_init, maxcv_init);

    % Determine the least merit value for each problem.
    merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
    if ismember(feature_name, {FeatureName.NOISY.value, FeatureName.TOUGH.value, FeatureName.TRUNCATED.value})
        feature_plain = Feature(FeatureName.PLAIN.value);
        fprintf("INFO: Starting the computation of the plain profiles.\n");
        [fun_histories_plain, maxcv_histories_plain] = solveAll(problem_names, problem_options, solvers, labels, feature_plain, max_eval_factor, profile_options);
        merit_histories_plain = computeMeritValues(fun_histories_plain, maxcv_histories_plain);
        merit_min_plain = min(min(min(merit_histories_plain, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
        merit_min = min(merit_min, merit_min_plain, 'omitnan');
    end

    % Paths to the results.
    timestamp = datestr(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH-mm-SSZ'), 'yyyy-mm-ddTHH-MM-SSZ');
    full_path = mfilename('fullpath');
    [folder_path, ~, ~] = fileparts(full_path);
    root_path = fileparts(folder_path);
    path_out = fullfile(root_path, 'out', feature.name, profile_options.(ProfileOptionKey.BENCHMARK_ID.value), timestamp);
    path_perf_out = fullfile(path_out, 'figures', 'perf');
    path_data_out = fullfile(path_out, 'figures', 'data');
    path_log_out = fullfile(path_out, 'figures', 'log');
    path_perf_hist_out = fullfile(path_perf_out, 'hist');
    path_data_hist_out = fullfile(path_data_out, 'hist');
    path_log_hist_out = fullfile(path_log_out, 'hist');
    path_perf_ret_out = fullfile(path_perf_out, 'ret');
    path_data_ret_out = fullfile(path_data_out, 'ret');
    path_log_ret_out = fullfile(path_log_out, 'ret');
    if ~exist(path_out, 'dir')
        mkdir(path_out);
    end
    if ~exist(path_perf_hist_out, 'dir')
        mkdir(path_perf_hist_out);
    end
    if ~exist(path_data_hist_out, 'dir')
        mkdir(path_data_hist_out);
    end
    if ~exist(path_perf_ret_out, 'dir')
        mkdir(path_perf_ret_out);
    end
    if ~exist(path_data_ret_out, 'dir')
        mkdir(path_data_ret_out);
    end


    % Store the names of the problems.
    path_txt = fullfile(path_out, 'problems.txt');
    fid = fopen(path_txt, 'w');
    if fid == -1
        error("Cannot open the file %s.", path_txt);
    end
    fprintf(fid, '%s\n', sort(string(problem_names)));
    fclose(fid);

    % Set the default values for plotting.
    fprintf("Creating results.\n");
    set(groot, 'DefaultLineLineWidth', 1);
    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');
    set(groot, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(groot, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    [n_problems, n_solvers, n_runs, max_eval] = size(merit_histories);
    tolerances = 10.^(-1:-1:-10);
    pdf_summary = 'summary.pdf';
    pdf_perf_hist = 'perf_hist.pdf';
    pdf_perf_ret = 'perf_ret.pdf';
    pdf_data_hist = 'data_hist.pdf';
    pdf_data_ret = 'data_ret.pdf';
    pdf_log_hist = 'log-ratio_hist.pdf';
    pdf_log_ret = 'log-ratio_ret.pdf';
    for i_profile = 1:10
        tolerance = tolerances(i_profile);
        fprintf("Creating profiles for tolerance %g.\n", tolerance);
        tolerance_label = ['$\tau = 10^{', int2str(log10(tolerance)), '}$'];

        work_hist = NaN(n_problems, n_solvers, n_runs);
        work_out = NaN(n_problems, n_solvers, n_runs);
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
                    if merit_out(i_problem, i_solver, i_run) <= threshold
                        work_out(i_problem, i_solver, i_run) = n_eval(i_problem, i_solver, i_run);
                    end
                end
            end
        end

        % Draw the profiles.
        










        % Calculate the x-axes of the performance and data profiles.
        [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data] = getExtendedPerformancesDataProfileAxes(work_hist, problem_dimensions);

        

        % Plot the performance profiles.
        fprintf("Creating performance profiles for tolerance %g.\n", tolerance);
        [fig, ~] = drawProfile(x_perf, y_perf, 'perf', ratio_max_perf, labels, tolerance_label);
        eps_filename_perf = fullfile(path_perf_hist_out, ['performance_profile_' int2str(i_profile) '.eps']);
        print(fig, eps_filename_perf, '-depsc');
        pdf_filename_perf = fullfile(path_perf_hist_out, ['performance_profile_' int2str(i_profile) '.pdf']);
        print(fig, pdf_filename_perf, '-dpdf');
        merge_pdf_filename_perf = fullfile(path_perf_hist_out, pdf_perf_hist);
        if i_profile == 1
            exportgraphics(fig, merge_pdf_filename_perf, 'ContentType', 'vector');
        else
            exportgraphics(fig, merge_pdf_filename_perf, 'ContentType', 'vector', 'Append', true);
        end
        close(fig);
        
        % Plot the data profiles.
        fprintf("Creating data profiles for tolerance %g.\n", tolerance);
        [fig, ~] = drawProfile(x_data, y_data, 'data', ratio_max_data, labels, tolerance_label);
        eps_filename_data = fullfile(path_data_hist_out, ['data_profile_' int2str(i_profile) '.eps']);
        print(fig, eps_filename_data, '-depsc');
        pdf_filename_data = fullfile(path_data_hist_out, ['data_profile_' int2str(i_profile) '.pdf']);
        print(fig, pdf_filename_data, '-dpdf');
        merge_pdf_filename_data = fullfile(path_data_hist_out, pdf_data_hist);
        if i_profile == 1
            exportgraphics(fig, merge_pdf_filename_data, 'ContentType', 'vector');
        else
            exportgraphics(fig, merge_pdf_filename_data, 'ContentType', 'vector', 'Append', true);
        end
        close(fig);

        % Draw the log-ratio profiles.
        if n_solvers == 2

            if ~exist(path_log_hist_out, 'dir')
                mkdir(path_log_hist_out);
            end
            if ~exist(path_log_ret_out, 'dir')
                mkdir(path_log_ret_out);
            end

            work_flat = reshape(permute(work_hist, [1, 3, 2]), n_problems * n_runs, n_solvers);
            log_ratio = NaN(n_problems * n_runs, 1);
            log_ratio_finite = isfinite(work_flat(:, 1)) & isfinite(work_flat(:, 2));
            log_ratio(log_ratio_finite) = log2(work_flat(log_ratio_finite, 1)) - log2(work_flat(log_ratio_finite, 2));
            ratio_max = max(max(abs(log_ratio(log_ratio_finite)), [], 'all'), eps);
            log_ratio(isnan(work_flat(:, 1)) & isfinite(work_flat(:, 2))) = 2.0 * ratio_max;
            log_ratio(isfinite(work_flat(:, 1)) & isnan(work_flat(:, 2))) = -2.0 * ratio_max;
            log_ratio(isnan(work_flat(:, 1)) & isnan(work_flat(:, 2))) = 0.0;
            log_ratio = sort(log_ratio);

            fig = figure('Visible', 'off');
            x = [1:(n_problems * n_runs)]';
            bar(x(log_ratio < 0), log_ratio(log_ratio < 0));
            hold on;
            bar(x(log_ratio > 0), log_ratio(log_ratio > 0));
            text((n_problems * n_runs + 1) / 2, -ratio_max, labels{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text((n_problems * n_runs + 1) / 2, ratio_max, labels{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            xticks([]);
            xlim([0.5, n_problems * n_runs + 0.5]);
            ylim([-1.1 * ratio_max, 1.1 * ratio_max]);
            xlabel('Problem', 'Interpreter', 'latex');
            ylabel(['Log-ratio profile ' tolerance_label], 'Interpreter', 'latex');
            eps_filename_log = fullfile(path_log_hist_out, ['log_ratio_profile_' int2str(i_profile) '.eps']);
            print(fig, eps_filename_log, '-depsc');
            pdf_filename_log = fullfile(path_log_hist_out, ['log_ratio_profile_' int2str(i_profile) '.pdf']);
            print(fig, pdf_filename_log, '-dpdf');
            merge_pdf_filename_log = fullfile(path_log_hist_out, pdf_log_hist);
            if i_profile == 1
                exportgraphics(fig, merge_pdf_filename_log, 'ContentType', 'vector');
            else
                exportgraphics(fig, merge_pdf_filename_log, 'ContentType', 'vector', 'Append', true);
            end
            close(fig);
        end
    end

    % Move the merged pdf files to the "path_out" folder.
    movefile(merge_pdf_filename_perf, fullfile(path_perf_out, pdf_perf_hist));
    movefile(merge_pdf_filename_data, fullfile(path_data_out, pdf_data_hist));
    if n_solvers == 2
        movefile(merge_pdf_filename_log, fullfile(path_log_out, pdf_log_hist));
    end

end


function merit_histories = computeMeritValues(fun_histories, maxcv_histories)
    is_nearly_feasible = maxcv_histories <= 1e-12;
    is_very_infeasible = maxcv_histories >= 1e-6;
    is_undecided = ~is_nearly_feasible & ~is_very_infeasible;
    merit_histories = NaN(size(fun_histories));
    merit_histories(is_nearly_feasible) = fun_histories(is_nearly_feasible);
    merit_histories(is_very_infeasible) = Inf;
    merit_histories(is_undecided) = fun_histories(is_undecided) + 1e8 * maxcv_histories(is_undecided);
end