function [results_plibs, profile_options] = loadResults(problem_options, profile_options)
%LOADRESULTS loads the results by the given options.

    results_plibs = {};
    if ~isfield(profile_options, ProfileOptionKey.LOAD.value) || isempty(profile_options.(ProfileOptionKey.LOAD.value))
        return;
    end

    % Check whether it is 'latest' or a string (or char) in the format of 'yyyyMMdd_HHmmss'.
    time_stamp_pattern = '^\d{4}(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})$';
    if ~strcmp(profile_options.(ProfileOptionKey.LOAD.value), 'latest') && isempty(regexp(profile_options.(ProfileOptionKey.LOAD.value), time_stamp_pattern, 'once'))
        error("MATLAB:checkValidityProfileOptions:loadNotValid", "The option `load` should be either 'latest' or a time stamp in the format of 'yyyyMMdd_HHmmss'.");
    end

    % Set the path to search for the data to load.
    if isfield(profile_options, ProfileOptionKey.BENCHMARK_ID.value)
        search_path = fullfile(pwd, profile_options.(ProfileOptionKey.BENCHMARK_ID.value));
    else
        search_path = pwd;
        profile_options.(ProfileOptionKey.BENCHMARK_ID.value) = '.';
    end

    % Find the path of the data to load.
    time_stamp_files = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
    if strcmp(profile_options.(ProfileOptionKey.LOAD.value), 'latest')
        % Try to find all the txt files named by time_stamp in the search_path directory and find the latest one.
        % Note that we limit the search to 5 levels of subdirectories.
        time_stamp_files = search_in_dir(search_path, 'time_stamp_*.txt', 5, 0, time_stamp_files);
        if isempty(time_stamp_files)
            error("MATLAB:loadResults:NoTimeStamps", "Failed to load data since no time_stamp files are found in the directory '%s'. Note that the search is limited to 5 levels of subdirectories.", search_path);
        end
        time_stamps = arrayfun(@(f) datetime(f.name(12:end-4), 'InputFormat', 'yyyyMMdd_HHmmss'), time_stamp_files);
        [~, indexes] = sort(time_stamps, 'descend');
        % Get the latest time_stamp file. If there are multiple files with the same time_stamp, we take the first one.
        latest_idx = indexes(1);
        latest_time_stamp_file = time_stamp_files(latest_idx);
        path_data = latest_time_stamp_file.folder;
    else
        % Same as above, but we only search for the specific time_stamp file.
        pattern = ['time_stamp_', profile_options.(ProfileOptionKey.LOAD.value), '.txt'];
        time_stamp_file = search_in_dir(search_path, pattern, 5, 0, time_stamp_files);
        if isempty(time_stamp_file)
            error("MATLAB:loadResults:NoTimeStamps", "Failed to load data since no time_stamp named '%s' is found in the directory '%s'. Note that the search is limited to 5 levels of subdirectories.", ['time_stamp_', profile_options.(ProfileOptionKey.LOAD.value), '.txt'], search_path);
        end
        path_data = time_stamp_file.folder;
    end

    % Load data from the 'data_for_loading.mat' file in the path_data directory.
    warning('off');
    fprintf("\nINFO: Loading data from the directory '%s'...\n", path_data);
    load(fullfile(path_data, 'data_for_loading.mat'), 'results_plibs');
    warning('on');
    if ~exist('results_plibs', 'var')
        error("MATLAB:loadResults:NoResultsPlibsDataMatFile", "Failed to load the variable 'results_plibs' from the 'data_for_loading.mat' file in the directory '%s'.", path_data);
    end

    % Process the 'solvers_to_load' field.
    n_solvers_loaded = size(results_plibs{1}.fun_histories, 2);
    if isfield(profile_options, ProfileOptionKey.SOLVERS_TO_LOAD.value)
        solvers_to_load = profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value);
        if any(solvers_to_load > n_solvers_loaded) || (isfield(profile_options, ProfileOptionKey.SOLVER_NAMES.value) && numel(profile_options.(ProfileOptionKey.SOLVER_NAMES.value)) ~= numel(solvers_to_load))
            error("MATLAB:loadResults:solvers_to_loadNotValid", "The option `solvers_to_load` should be a vector of different integers selected from 1 to the total number of solvers in the loaded data, and at least two indices should be provided.");
        end
    else
        solvers_to_load = 1:n_solvers_loaded;
    end

    % Assign values to 'profile_options' by the loaded 'results_plibs'.
    if ~isfield(profile_options, ProfileOptionKey.SOLVER_NAMES.value)
        profile_options.(ProfileOptionKey.SOLVER_NAMES.value) = results_plibs{1}.solver_names(solvers_to_load);
    end
    if ~isfield(profile_options, ProfileOptionKey.FEATURE_STAMP.value)
        profile_options.(ProfileOptionKey.FEATURE_STAMP.value) = results_plibs{1}.feature_stamp;
    end

    % Select problem libraries to load by the 'plibs' field.
    if isfield(problem_options, ProblemOptionKey.PLIBS.value)
        valid_plib_names = cellfun(@(x) x.plib, results_plibs, 'UniformOutput', false);
        plibs = problem_options.(ProblemOptionKey.PLIBS.value);
        if ~all(cellfun(@ischarstr, plibs)) || ~all(ismember(plibs, valid_plib_names))
            error("MATLAB:loadResults:plibNotValid", "When you use the option `load`, the option `plibs` should be a cell array of strings or chars, and each string or char should be one of the problem libraries in the loaded data ('%s').", strjoin(valid_plib_names, "', '"));
        end
        % Select the problem libraries to load.
        results_plibs = results_plibs(cellfun(@(x) ismember(x.plib, plibs), results_plibs));
        if isempty(results_plibs)
            error("MATLAB:loadResults:NoValidPlibs", "No problem libraries are selected from the loaded data.");
        end
    end

    % Truncate data by the 'solvers_to_load' field.
    for i_plib = 1:size(results_plibs, 2)
        results_plib = results_plibs{i_plib};
        results_plib = truncate_solvers(results_plib, solvers_to_load);
        if isfield(results_plib, 'results_plib_plain')
            results_plib.results_plib_plain = truncate_solvers(results_plib.results_plib_plain, solvers_to_load);
        end
        results_plibs{i_plib} = results_plib;
    end

    % Truncate data by 'problem_options'.
    for i_plib = 1:size(results_plibs, 2)
        results_plib = results_plibs{i_plib};
        results_plib = truncate_problems(results_plib, problem_options);
        if isfield(results_plib, 'results_plib_plain')
            results_plib.results_plib_plain = truncate_problems(results_plib.results_plib_plain, problem_options);
        end
        results_plibs{i_plib} = results_plib;
    end

    % Recompute merit values if profile_options.merit_fun is provided.
    if isfield(profile_options, ProfileOptionKey.MERIT_FUN.value)
        merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
        for i_plib = 1:size(results_plibs, 2)
            results_plib = results_plibs{i_plib};
            results_plib = recompute_merits(results_plib, merit_fun);
            if isfield(results_plib, 'results_plib_plain')
                results_plib.results_plib_plain = recompute_merits(results_plib.results_plib_plain, merit_fun);
            end
            results_plibs{i_plib} = results_plib;
        end
    end
end

function files = search_in_dir(current_path, pattern, max_depth, current_depth, files)
    % Recursive function to search for files with a specified pattern within a limited directory depth.

    % Check if the current depth exceeds the maximum depth
    if current_depth > max_depth
        return; % Stop searching deeper
    end

    % Get files in the current directory matching the pattern
    file_list = dir(fullfile(current_path, pattern));
    for i = 1:length(file_list)
        if ~file_list(i).isdir
            % Append the file information to the result list
            files = [files; file_list(i)];
        end
    end

    % Get all subdirectories in the current directory
    sub_dirs = dir(current_path);
    for i = 1:length(sub_dirs)
        % Skip '.' and '..' directories
        if sub_dirs(i).isdir && ~strcmp(sub_dirs(i).name, '.') && ~strcmp(sub_dirs(i).name, '..')
            % Recursively search in the subdirectory
            files = search_in_dir(fullfile(current_path, sub_dirs(i).name), pattern, max_depth, current_depth + 1, files);
        end
    end
end

function results_plib = truncate_solvers(results_plib, solvers_to_load)
    % Truncate the loaded data by the 'solvers_to_load' field.
    results_plib.solver_names = results_plib.solver_names(solvers_to_load);
    results_plib.fun_histories = results_plib.fun_histories(:, solvers_to_load, :, :);
    results_plib.maxcv_histories = results_plib.maxcv_histories(:, solvers_to_load, :, :);
    results_plib.fun_outs = results_plib.fun_outs(:, solvers_to_load, :);
    results_plib.maxcv_outs = results_plib.maxcv_outs(:, solvers_to_load, :);
    results_plib.n_evals = results_plib.n_evals(:, solvers_to_load, :);
    results_plib.computation_times = results_plib.computation_times(:, solvers_to_load, :);
    results_plib.solvers_successes = results_plib.solvers_successes(:, solvers_to_load, :);
    results_plib.merit_histories = results_plib.merit_histories(:, solvers_to_load, :, :);
    results_plib.merit_outs = results_plib.merit_outs(:, solvers_to_load, :);
end

function results_plib = truncate_problems(results_plib, problem_options)
    % Truncate the loaded data by the 'problem_options' field.
    p_to_load = true(size(results_plib.problem_names));

    if isfield(problem_options, ProblemOptionKey.PTYPE.value)
        % Judge whether the type of the problems in the loaded data satisfies the options.
        results_plib.ptype = intersect(results_plib.ptype, problem_options.(ProblemOptionKey.PTYPE.value));
        for i = 1:numel(results_plib.problem_types)
            if ~ismember(results_plib.problem_types{i}, problem_options.(ProblemOptionKey.PTYPE.value))
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MINDIM.value)
        % Judge whether the dimension of the problems in the loaded data satisfies the options.
        results_plib.mindim = max(results_plib.mindim, problem_options.(ProblemOptionKey.MINDIM.value));
        for i = 1:numel(results_plib.problem_dims)
            if results_plib.problem_dims(i) < problem_options.(ProblemOptionKey.MINDIM.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MAXDIM.value)
        % Judge whether the dimension of the problems in the loaded data satisfies the options.
        results_plib.maxdim = min(results_plib.maxdim, problem_options.(ProblemOptionKey.MAXDIM.value));
        for i = 1:numel(results_plib.problem_dims)
            if results_plib.problem_dims(i) > problem_options.(ProblemOptionKey.MAXDIM.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MINB.value)
        % Judge whether the number of bound constraints of the problems in the loaded data
        % satisfies the options.
        results_plib.minb = max(results_plib.minb, problem_options.(ProblemOptionKey.MINB.value));
        for i = 1:numel(results_plib.problem_mbs)
            if results_plib.problem_mbs(i) < problem_options.(ProblemOptionKey.MINB.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MAXB.value)
        % Judge whether the number of bound constraints of the problems in the loaded data
        % satisfies the options.
        results_plib.maxb = min(results_plib.maxb, problem_options.(ProblemOptionKey.MAXB.value));
        for i = 1:numel(results_plib.problem_mbs)
            if results_plib.problem_mbs(i) > problem_options.(ProblemOptionKey.MAXB.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MINLCON.value)
        % Judge whether the number of linear constraints of the problems in the loaded data
        % satisfies the options.
        results_plib.minlcon = max(results_plib.minlcon, problem_options.(ProblemOptionKey.MINLCON.value));
        for i = 1:numel(results_plib.problem_mlcons)
            if results_plib.problem_mlcons(i) < problem_options.(ProblemOptionKey.MINLCON.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MAXLCON.value)
        % Judge whether the number of linear constraints of the problems in the loaded data
        % satisfies the options.
        results_plib.maxlcon = min(results_plib.maxlcon, problem_options.(ProblemOptionKey.MAXLCON.value));
        for i = 1:numel(results_plib.problem_mlcons)
            if results_plib.problem_mlcons(i) > problem_options.(ProblemOptionKey.MAXLCON.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MINNLCON.value)
        % Judge whether the number of nonlinear constraints of the problems in the loaded data
        % satisfies the options.
        results_plib.minnlcon = max(results_plib.minnlcon, problem_options.(ProblemOptionKey.MINNLCON.value));
        for i = 1:numel(results_plib.problem_mnlcons)
            if results_plib.problem_mnlcons(i) < problem_options.(ProblemOptionKey.MINNLCON.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MAXNLCON.value)
        % Judge whether the number of nonlinear constraints of the problems in the loaded data
        % satisfies the options.
        results_plib.maxnlcon = min(results_plib.maxnlcon, problem_options.(ProblemOptionKey.MAXNLCON.value));
        for i = 1:numel(results_plib.problem_mnlcons)
            if results_plib.problem_mnlcons(i) > problem_options.(ProblemOptionKey.MAXNLCON.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MINCON.value)
        % Judge whether the number of constraints of the problems in the loaded data satisfies
        % the options.
        results_plib.mincon = max(results_plib.mincon, problem_options.(ProblemOptionKey.MINCON.value));
        for i = 1:numel(results_plib.problem_mcons)
            if results_plib.problem_mcons(i) < problem_options.(ProblemOptionKey.MINCON.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.MAXCON.value)
        % Judge whether the number of constraints of the problems in the loaded data satisfies
        % the options.
        results_plib.maxcon = min(results_plib.maxcon, problem_options.(ProblemOptionKey.MAXCON.value));
        for i = 1:numel(results_plib.problem_mcons)
            if results_plib.problem_mcons(i) > problem_options.(ProblemOptionKey.MAXCON.value)
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.PROBLEM_NAMES.value)
        % Only load the problems whose names are in the 'problem_names' field.
        for i = 1:numel(results_plib.problem_names)
            if ~ismember(results_plib.problem_names{i}, problem_options.(ProblemOptionKey.PROBLEM_NAMES.value))
                p_to_load(i) = false;
            end
        end
    end
    if isfield(problem_options, ProblemOptionKey.EXCLUDELIST.value)
        % Judge whether the problems in the loaded data are in the exclude list.
        for i = 1:numel(results_plib.problem_names)
            if ismember(results_plib.problem_names{i}, problem_options.(ProblemOptionKey.EXCLUDELIST.value))
                p_to_load(i) = false;
            end
        end
    end

    % Check whether there is any problem to load.
    if ~any(p_to_load)
        fprintf("\nINFO: No problem is selected to load by the given options. Please check your options and try again.\n");
        results_plib = [];
        return;
    end

    % Truncate the loaded data.
    results_plib.fun_histories = results_plib.fun_histories(p_to_load, :, :, :);
    results_plib.maxcv_histories = results_plib.maxcv_histories(p_to_load, :, :, :);
    results_plib.fun_outs = results_plib.fun_outs(p_to_load, :, :);
    results_plib.maxcv_outs = results_plib.maxcv_outs(p_to_load, :, :);
    results_plib.fun_inits = results_plib.fun_inits(p_to_load, :);
    results_plib.maxcv_inits = results_plib.maxcv_inits(p_to_load, :);
    results_plib.n_evals = results_plib.n_evals(p_to_load, :, :);
    results_plib.problem_names = results_plib.problem_names(p_to_load);
    results_plib.problem_types = results_plib.problem_types(p_to_load);
    results_plib.problem_dims = results_plib.problem_dims(p_to_load);
    results_plib.problem_mbs = results_plib.problem_mbs(p_to_load);
    results_plib.problem_mlcons = results_plib.problem_mlcons(p_to_load);
    results_plib.problem_mnlcons = results_plib.problem_mnlcons(p_to_load);
    results_plib.problem_mcons = results_plib.problem_mcons(p_to_load);
    results_plib.computation_times = results_plib.computation_times(p_to_load, :, :);
    results_plib.solvers_successes = results_plib.solvers_successes(p_to_load, :, :);
    results_plib.merit_histories = results_plib.merit_histories(p_to_load, :, :, :);
    results_plib.merit_outs = results_plib.merit_outs(p_to_load, :, :);
    results_plib.merit_inits = results_plib.merit_inits(p_to_load, :);
end

function results_plib = recompute_merits(results_plib, merit_fun)
    % Recompute the merit values for the loaded data if the 'merit_fun' is provided.
    
    try
        results_plib.merit_histories = meritFunCompute(merit_fun, results_plib.fun_histories, results_plib.maxcv_histories, results_plib.maxcv_inits);
        results_plib.merit_outs = meritFunCompute(merit_fun, results_plib.fun_outs, results_plib.maxcv_outs, results_plib.maxcv_inits);
        results_plib.merit_inits = meritFunCompute(merit_fun, results_plib.fun_inits, results_plib.maxcv_inits, results_plib.maxcv_inits);
    catch
        error("MATLAB:loadResults:meritFunError", "The merit function provided in the options is not valid. Please check the function and try again.");
    end
end