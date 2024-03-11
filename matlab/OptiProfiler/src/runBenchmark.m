function runBenchmark(solvers, labels, problem_names, feature_name, varargin)
%RUNBENCHMARK creates the benchmark profiles.
%
%   The benchmark profiles include the performance and data profiles [1]_,
%   [2]_, [4]_, and the log-ratio profiles [3]_, [5]_. The log-ratio profiles
%   are available only when there are exactly two solvers.
%
%   References:
%   .. [1] E. D. Dolan and J. J. Moré. Benchmarking optimization software with
%          performance profiles. *Math. Program.*, 91(2):201–213, 2002.
%          `doi:10.1007/s101070100263
%          <https://doi.org/10.1007/s101070100263>.
%   .. [2] N. Gould and J. Scott. A note on performance profiles for
%          benchmarking software. *ACM Trans. Math. Software*, 43(2):15:1–5,
%          2016. `doi:10.1145/2950048 <https://doi.org/10.1145/2950048>.
%   .. [3] J. L. Morales. A numerical study of limited memory BFGS methods.
%          *Appl. Math. Lett.*, 15(4):481–487, 2002.
%          `doi:10.1016/S0893-9659(01)00162-8
%          <https://doi.org/10.1016/S0893-9659(01)00162-8>.
%   .. [4] J. J. Moré and S. M. Wild. Benchmarking derivative-free optimization
%          algorithms. *SIAM J. Optim.*, 20(1):172–191, 2009.
%          `doi:10.1137/080724083 <https://doi.org/10.1137/080724083>.
%   .. [5] H.-J. M. Shi, M. Q. Xuan, F. Oztoprak, and J. Nocedal. On the
%          numerical performance of finite-difference-based methods for
%          derivative-free optimization. *Optim. Methods Softw.*,
%          38(2):289–311, 2023. `doi:10.1080/10556788.2022.2121832
%          <https://doi.org/10.1080/10556788.2022.2121832>.


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
    [fun_histories, maxcv_histories, fun_ret, maxcv_ret, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, time_processes] = solveAllProblems(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);
    merit_histories = computeMeritValues(fun_histories, maxcv_histories, maxcv_init);
    merit_ret = computeMeritValues(fun_ret, maxcv_ret, maxcv_init);
    merit_init = computeMeritValues(fun_init, maxcv_init, maxcv_init);

    % Determine the least merit value for each problem.
    merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
    if feature.isStochastic
        feature_plain = Feature(FeatureName.PLAIN.value);
        fprintf("INFO: Starting the computation of the plain profiles.\n");
        [fun_histories_plain, maxcv_histories_plain, ~, ~, ~, ~, ~, ~, ~, time_processes_plain] = solveAllProblems(problem_names, problem_options, solvers, labels, feature_plain, max_eval_factor, profile_options);
        time_processes = time_processes + time_processes_plain;
        merit_histories_plain = computeMeritValues(fun_histories_plain, maxcv_histories_plain, maxcv_init);
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
    path_log_ratio_out = fullfile(path_out, 'figures', 'log-ratio');
    path_perf_hist_out = fullfile(path_perf_out, 'hist');
    path_data_hist_out = fullfile(path_data_out, 'hist');
    path_log_ratio_hist_out = fullfile(path_log_ratio_out, 'hist');
    path_perf_ret_out = fullfile(path_perf_out, 'ret');
    path_data_ret_out = fullfile(path_data_out, 'ret');
    path_log_ratio_ret_out = fullfile(path_log_ratio_out, 'ret');
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
    [sorted_problem_names, idx] = sort(problem_names);
    sorted_time_processes = time_processes(idx);
    fid = fopen(path_txt, 'w');
    if fid == -1
        error("Cannot open the file %s.", path_txt);
    end
    for i = 1:length(sorted_problem_names)
        fprintf(fid, "%s: %.2f seconds\n", sorted_problem_names{i}, sorted_time_processes(i));
    end
    fclose(fid);

    % Set the default values for plotting.
    fprintf("Creating results.\n");
    set(groot, 'DefaultLineLineWidth', 1);
    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');
    set(groot, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(groot, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    [n_problems, n_solvers, n_runs, max_eval] = size(merit_histories);
    if n_solvers <= 2
        if ~exist(path_log_ratio_hist_out, 'dir')
            mkdir(path_log_ratio_hist_out);
        end
        if ~exist(path_log_ratio_ret_out, 'dir')
            mkdir(path_log_ratio_ret_out);
        end
    end
    tolerances = 10.^(-1:-1:-10);
    pdf_summary = fullfile(path_out, 'summary.pdf');
    pdf_perf_hist_summary = fullfile(path_perf_out, 'perf_hist.pdf');
    pdf_perf_ret_summary = fullfile(path_perf_out, 'perf_ret.pdf');
    pdf_data_hist_summary = fullfile(path_data_out, 'data_hist.pdf');
    pdf_data_ret_summary = fullfile(path_data_out, 'data_ret.pdf');
    pdf_log_ratio_hist_summary = fullfile(path_log_ratio_out, 'log-ratio_hist.pdf');
    pdf_log_ratio_ret_summary = fullfile(path_log_ratio_out, 'log-ratio_ret.pdf');
    for i_profile = 1:10
        tolerance = tolerances(i_profile);
        fprintf("Creating profiles for tolerance %g.\n", tolerance);
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

        % Draw the profiles.
        [fig_summary, fig_perf_hist, fig_perf_ret, fig_data_hist, fig_data_ret, fig_log_ratio_hist, fig_log_ratio_ret] = drawProfiles(work_hist, work_ret, problem_dimensions, labels, tolerance_label);
        eps_perf_hist = fullfile(path_perf_hist_out, ['perf_hist_' int2str(i_profile) '.eps']);
        print(fig_perf_hist, eps_perf_hist, '-depsc');
        pdf_perf_hist = fullfile(path_perf_hist_out, ['perf_hist_' int2str(i_profile) '.pdf']);
        print(fig_perf_hist, pdf_perf_hist, '-dpdf');
        eps_perf_ret = fullfile(path_perf_ret_out, ['perf_ret_' int2str(i_profile) '.eps']);
        print(fig_perf_ret, eps_perf_ret, '-depsc');
        pdf_perf_ret = fullfile(path_perf_ret_out, ['perf_ret_' int2str(i_profile) '.pdf']);
        print(fig_perf_ret, pdf_perf_ret, '-dpdf');
        eps_data_hist = fullfile(path_data_hist_out, ['data_hist_' int2str(i_profile) '.eps']);
        print(fig_data_hist, eps_data_hist, '-depsc');
        pdf_data_hist = fullfile(path_data_hist_out, ['data_hist_' int2str(i_profile) '.pdf']);
        print(fig_data_hist, pdf_data_hist, '-dpdf');
        eps_data_ret = fullfile(path_data_ret_out, ['data_ret_' int2str(i_profile) '.eps']);
        print(fig_data_ret, eps_data_ret, '-depsc');
        pdf_data_ret = fullfile(path_data_ret_out, ['data_ret_' int2str(i_profile) '.pdf']);
        print(fig_data_ret, pdf_data_ret, '-dpdf');
        if ~isempty(fig_log_ratio_hist)
            eps_log_ratio_hist = fullfile(path_log_ratio_hist_out, ['log-ratio_hist_' int2str(i_profile) '.eps']);
            print(fig_log_ratio_hist, eps_log_ratio_hist, '-depsc');
            pdf_log_ratio_hist = fullfile(path_log_ratio_hist_out, ['log-ratio_hist_' int2str(i_profile) '.pdf']);
            print(fig_log_ratio_hist, pdf_log_ratio_hist, '-dpdf');
        end
        if ~isempty(fig_log_ratio_ret)
            eps_log_ratio_ret = fullfile(path_log_ratio_ret_out, ['log-ratio_ret_' int2str(i_profile) '.eps']);
            print(fig_log_ratio_ret, eps_log_ratio_ret, '-depsc');
            pdf_log_ratio_ret = fullfile(path_log_ratio_ret_out, ['log-ratio_ret_' int2str(i_profile) '.pdf']);
            print(fig_log_ratio_ret, pdf_log_ratio_ret, '-dpdf');
        end
        if i_profile == 1
            exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector');
            exportgraphics(fig_perf_hist, pdf_perf_hist_summary, 'ContentType', 'vector');
            exportgraphics(fig_perf_ret, pdf_perf_ret_summary, 'ContentType', 'vector');
            exportgraphics(fig_data_hist, pdf_data_hist_summary, 'ContentType', 'vector');
            exportgraphics(fig_data_ret, pdf_data_ret_summary, 'ContentType', 'vector');
            if ~isempty(fig_log_ratio_hist)
                exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist_summary, 'ContentType', 'vector');
            end
            if ~isempty(fig_log_ratio_ret)
                exportgraphics(fig_log_ratio_ret, pdf_log_ratio_ret_summary, 'ContentType', 'vector');
            end
        else
            exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector', 'Append', true);
            exportgraphics(fig_perf_hist, pdf_perf_hist_summary, 'ContentType', 'vector', 'Append', true);
            exportgraphics(fig_perf_ret, pdf_perf_ret_summary, 'ContentType', 'vector', 'Append', true);
            exportgraphics(fig_data_hist, pdf_data_hist_summary, 'ContentType', 'vector', 'Append', true);
            exportgraphics(fig_data_ret, pdf_data_ret_summary, 'ContentType', 'vector', 'Append', true);
            if ~isempty(fig_log_ratio_hist)
                exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist_summary, 'ContentType', 'vector', 'Append', true);
            end
            if ~isempty(fig_log_ratio_ret)
                exportgraphics(fig_log_ratio_ret, pdf_log_ratio_ret_summary, 'ContentType', 'vector', 'Append', true);
            end
        end

        % Close the figures.
        close(fig_summary);
        close(fig_perf_hist);
        close(fig_perf_ret);
        close(fig_data_hist);
        close(fig_data_ret);
        if ~isempty(fig_log_ratio_hist)
            close(fig_log_ratio_hist);
        end
        if ~isempty(fig_log_ratio_ret)
            close(fig_log_ratio_ret);
        end

    end

end


function merit_values = computeMeritValues(fun_values, maxcv_values, maxcv_init)
    copied_dim = size(fun_values);
    maxcv_init = repmat(maxcv_init, [1, copied_dim(2:end)]);
    infeasibility_thresholds = max(1e-5, maxcv_init);
    is_infeasible = maxcv_values > infeasibility_thresholds;
    is_almost_feasible = (1e-10 < maxcv_values) & (maxcv_values <= infeasibility_thresholds);
    merit_values = fun_values;
    merit_values(is_infeasible | isnan(merit_values)) = Inf;
    merit_values(is_almost_feasible) = merit_values(is_almost_feasible) + 1e5 * maxcv_values(is_almost_feasible);
end