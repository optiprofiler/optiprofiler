function createProfiles(solvers, labels, problem_names, feature_name, varargin)
    %CREATEPROFILES plots performance profiles and data profiles.

    % TODO: the current script version does not support options, only supports "simple" features.

    % Set default values for optional arguments.
    myCluster = parcluster('local');
    nb_cores = myCluster.NumWorkers;
    profile_options = struct(ProfileOptionKey.N_JOBS.value, nb_cores);
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

    % Build feature.
    % TODO: now we first use "simple" default feature.
    feature = Feature(feature_name);
    fprintf("INFO: Starting the computation of the %s profiles.\n", feature.name);

    % Solve all the problems.
    max_eval_factor = 500;
    [fun_values, maxcv_values, fun_inits, maxcv_inits, n_evals, problem_names, problem_dimensions] = solveAll(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);

    % Compute the merit values.
    merit_values = computeMeritValues(fun_values, maxcv_values);
    merit_inits = computeMeritValues(fun_inits, maxcv_inits);

    % Determine the least merit value for each problem.
    merit_min = min(min(min(merit_values, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
    if ismember(feature_name, {FeatureName.NOISY.value, FeatureName.TOUGH.value, FeatureName.TRUNCATED.value})
        feature_plain = Feature(FeatureName.PLAIN.value);
        fprintf("INFO: Starting the computation of the plain profiles.\n");
        [fun_values_plain, maxcv_values_plain] = solveAll(problem_names, problem_options, solvers, labels, feature_plain, max_eval_factor, profile_options);
        merit_values_plain = computeMeritValues(fun_values_plain, maxcv_values_plain);
        merit_min_plain = min(min(min(merit_values_plain, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
        merit_min = min(merit_min, merit_min_plain, 'omitnan');
    end

    % Paths to the results.
    full_path = mfilename('fullpath');
    [folder_path, ~, ~] = fileparts(full_path);
    root_path = fileparts(folder_path);
    path_out = fullfile(root_path, 'out', feature.name);
    path_perf_out = fullfile(path_out, 'figures', 'perf');
    path_data_out = fullfile(path_out, 'figures', 'data');
    if ~exist(path_out, 'dir')
        mkdir(path_out);
    end
    if ~exist(path_perf_out, 'dir')
        mkdir(path_perf_out);
    end
    if ~exist(path_data_out, 'dir')
        mkdir(path_data_out);
    end
    timestamp = datestr(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH-mm-SSZ'), 'yyyy-mm-ddTHH-MM-SSZ');
    name_pdf_merged_perf = ['performance_profiles_' timestamp '.pdf'];
    name_pdf_merged_data = ['data_profiles_' timestamp '.pdf'];
    



    % Create the performance and data profiles.
    fprintf("Creating the results.\n");
    set(groot, 'DefaultLineLineWidth', 1);
    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');
    set(groot, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(groot, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    [n_problems, n_solvers, n_runs, max_eval] = size(merit_values);
    tolerances = 10.^(-1:-1:-10);
    for i_profile = 1:10
        tolerance = tolerances(i_profile);

        work = NaN(n_problems, n_solvers, n_runs);
        for i_problem = 1:n_problems
            if isfinite(merit_min(i_problem))
                threshold = max(tolerance * merit_inits(i_problem) + (1 - tolerance) * merit_min(i_problem), merit_min(i_problem));
            else
                threshold = -Inf;
            end
            for i_solver = 1:n_solvers
                for i_run = 1:n_runs
                    if min(merit_values(i_problem, i_solver, i_run, :), [], 'omitnan') <= threshold
                        work(i_problem, i_solver, i_run) = find(merit_values(i_problem, i_solver, i_run, :) <= threshold, 1, 'first');
                    end
                end
            end
        end

        % Calculate the x-axes of the performance and data profiles.
        denominator_perf = @(i_problem, i_run) min(work(i_problem, :, i_run), [], 'omitnan');
        [x_perf, y_perf, ratio_max_perf] = profileAxes(work, denominator_perf);
        x_perf(isinf(x_perf)) = ratio_max_perf ^ 2;
        denominator_data = @(i_problem, i_run) problem_dimensions(i_problem) + 1;
        [x_data, y_data, ratio_max_data] = profileAxes(work, denominator_data);
        x_data(isinf(x_data)) = 2 * ratio_max_data;

        tolerance_label = ['$\tau = 10^{', int2str(log10(tolerance)), '}$'];

        % Plot the performance profiles.
        fprintf("Creating performance profiles for tolerance %g.\n", tolerance);
        [fig, ~] = drawProfile(x_perf, y_perf, 'perf', ratio_max_perf, labels, tolerance_label);
        eps_filename_perf = fullfile(path_perf_out, ['performance_profile_' int2str(i_profile) '.eps']);
        print(fig, eps_filename_perf, '-depsc');
        pdf_filename_perf = fullfile(path_perf_out, ['performance_profile_' int2str(i_profile) '.pdf']);
        print(fig, pdf_filename_perf, '-dpdf');
        merge_pdf_filename_perf = fullfile(path_perf_out, name_pdf_merged_perf);
        if i_profile == 1
            exportgraphics(fig, merge_pdf_filename_perf, 'ContentType', 'vector');
        else
            exportgraphics(fig, merge_pdf_filename_perf, 'ContentType', 'vector', 'Append', true);
        end
        close(fig);
        
        % Plot the data profiles.
        fprintf("Creating data profiles for tolerance %g.\n", tolerance);
        [fig, ~] = drawProfile(x_data, y_data, 'data', ratio_max_data, labels, tolerance_label);
        eps_filename_data = fullfile(path_data_out, ['data_profile_' int2str(i_profile) '.eps']);
        print(fig, eps_filename_data, '-depsc');
        pdf_filename_data = fullfile(path_data_out, ['data_profile_' int2str(i_profile) '.pdf']);
        print(fig, pdf_filename_data, '-dpdf');
        merge_pdf_filename_data = fullfile(path_data_out, name_pdf_merged_data);
        if i_profile == 1
            exportgraphics(fig, merge_pdf_filename_data, 'ContentType', 'vector');
        else
            exportgraphics(fig, merge_pdf_filename_data, 'ContentType', 'vector', 'Append', true);
        end
        close(fig);

    end

    % Move the merged pdf files to the "path_out" folder.
    movefile(merge_pdf_filename_perf, fullfile(path_out, name_pdf_merged_perf));
    movefile(merge_pdf_filename_data, fullfile(path_out, name_pdf_merged_data));
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