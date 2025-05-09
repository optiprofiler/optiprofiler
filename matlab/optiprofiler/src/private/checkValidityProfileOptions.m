function profile_options = checkValidityProfileOptions(solvers, profile_options)
%CHECKVALIDITYPROFILEOPTIONS Check the validity of the options in profile_options

    if exist('parcluster', 'file') == 2 % Check if Parallel Computing Toolbox is available
        myCluster = parcluster('local');
        % Get the number of workers in the cluster
        nb_cores = myCluster.NumWorkers;
    else
        nb_cores = 1;
    end


    % Judge whether profile_options.n_jobs is a integer between 1 and nb_cores.
    if isfield(profile_options, ProfileOptionKey.N_JOBS.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.N_JOBS.value))
            error("MATLAB:checkValidityProfileOptions:n_jobsNotValid", "The field 'n_jobs' of options should be a integer.");
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) < 1
            profile_options.(ProfileOptionKey.N_JOBS.value) = 1;
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) > nb_cores
            profile_options.(ProfileOptionKey.N_JOBS.value) = nb_cores;
        else
            profile_options.(ProfileOptionKey.N_JOBS.value) = round(profile_options.(ProfileOptionKey.N_JOBS.value));
        end
    end


    % Judge whether profile_options.keep_pool is a boolean.
    if isfield(profile_options, ProfileOptionKey.KEEP_POOL.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.KEEP_POOL.value))
            error("MATLAB:checkValidityProfileOptions:keep_poolNotValid", "The field 'keep_pool' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.seed is a positive integer.
    if isfield(profile_options, ProfileOptionKey.SEED.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.SEED.value)) || profile_options.(ProfileOptionKey.SEED.value) <= 0
            error("MATLAB:checkValidityProfileOptions:seedNotValid", "The field 'seed' of options should be a positive integer.");
        end
    end


    % Judge whether profile_options.benchmark_id is a char or a string and satisfies the file name requirements.
    if isfield(profile_options, ProfileOptionKey.BENCHMARK_ID.value)
        is_valid_foldername = @(x) ischarstr(x) && ~isempty(x) && all(ismember(char(x), ['a':'z', 'A':'Z', '0':'9', '_', '-', '.']));
        if ~ischarstr(profile_options.(ProfileOptionKey.BENCHMARK_ID.value)) || ~is_valid_foldername(profile_options.(ProfileOptionKey.BENCHMARK_ID.value))
            error("MATLAB:checkValidityProfileOptions:benchmark_idNotValid", "The field 'benchmark_id' of options should be a char or a string satisfying the strict file name requirements (only containing letters, numbers, underscores, hyphens, and dots).");
        end
    end


    % Judge whether profile_options.solver_names is a cell of chars or strings.
    if isfield(profile_options, ProfileOptionKey.SOLVER_NAMES.value)
        if ~iscell(profile_options.(ProfileOptionKey.SOLVER_NAMES.value)) || ~all(cellfun(@(l) ischarstr(l), profile_options.(ProfileOptionKey.SOLVER_NAMES.value)))
            error("MATLAB:checkValidityProfileOptions:solver_namesNotCellOfcharstr", "The field 'solver_names' of options must be a cell of chars or strings.");
        end
        if ~isempty(solvers) && numel(profile_options.(ProfileOptionKey.SOLVER_NAMES.value)) ~= 0 && numel(profile_options.(ProfileOptionKey.SOLVER_NAMES.value)) ~= numel(solvers)
            error("MATLAB:checkValidityProfileOptions:solver_namesAndsolversLengthNotSame", "The number of the field 'solver_names' of options must equal the number of solvers.");
        end
        if numel(profile_options.(ProfileOptionKey.SOLVER_NAMES.value)) == 0 && ~isempty(solvers)
            profile_options.(ProfileOptionKey.SOLVER_NAMES.value) = cellfun(@(s) func2str(s), solvers, 'UniformOutput', false);
        end
        % Handle the case where the solver names are not valid MATLAB variable names.
        % Replace underscores with backslashes.
        profile_options.(ProfileOptionKey.SOLVER_NAMES.value) = cellfun(@(s) strrep(s, '_', '\_'), profile_options.(ProfileOptionKey.SOLVER_NAMES.value), 'UniformOutput', false);
    end


    % Judge whether profile_options.solver_isrand is a logical array of the same length as the number of solvers.
    if isfield(profile_options, ProfileOptionKey.SOLVER_ISRAND.value)
        if ~islogical(profile_options.(ProfileOptionKey.SOLVER_ISRAND.value)) || (~isempty(solvers) && numel(profile_options.(ProfileOptionKey.SOLVER_ISRAND.value)) ~= numel(solvers))
            error("MATLAB:checkValidityProfileOptions:solver_israndNotLogical", "The field 'solver_isrand' of options must be a logical array of the same length as the number of solvers.");
        end
    end


    % Judge whether profile_options.feature_stamp is a char or a string and satisfies the file name requirements.
    if isfield(profile_options, ProfileOptionKey.FEATURE_STAMP.value)
        if ~ischarstr(profile_options.(ProfileOptionKey.FEATURE_STAMP.value))
            error("MATLAB:checkValidityProfileOptions:feature_stampNotcharstr", "The field 'feature_stamp' of options must be a char or string.");
        end
        is_valid_foldername = @(x) ischarstr(x) && ~isempty(x) && all(ismember(char(x), ['a':'z', 'A':'Z', '0':'9', '_', '-', '.']));
        if ~is_valid_foldername(profile_options.(ProfileOptionKey.FEATURE_STAMP.value))
            error("MATLAB:checkValidityProfileOptions:feature_stampNotValid", "The field 'feature_stamp' of options should be a char or a string satisfying the strict file name requirements (only containing letters, numbers, underscores, hyphens, and dots).");
        end
    end


    % Judge whether profile_options.errorbar_type is among 'minmax' and 'meanstd'.
    if isfield(profile_options, ProfileOptionKey.ERRORBAR_TYPE.value)
        if ~ischarstr(profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value)) || ~ismember(char(profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value)), {'minmax', 'meanstd'})
            error("MATLAB:checkValidityProfileOptions:errorbar_typeNotValid", "The field 'errorbar_type' of options should be either 'minmax' or 'meanstd'.");
        end
    end


    % Judge whether profile_options.savepath is a string and exists. If not exists, create it.
    if isfield(profile_options, ProfileOptionKey.SAVEPATH.value)
        if ~ischarstr(profile_options.(ProfileOptionKey.SAVEPATH.value))
            error("MATLAB:checkValidityProfileOptions:savepathNotValid", "The field 'savepath' of options should be a char or a string.");
        elseif ~exist(profile_options.(ProfileOptionKey.SAVEPATH.value), 'dir')
            status = mkdir(profile_options.(ProfileOptionKey.SAVEPATH.value));
            if ~status
                error("MATLAB:checkValidityProfileOptions:savepathNotExist", "The field 'savepath' of options does not exist and cannot be created.");
            end
        end
    end


    % Judge whether profile_options.max_tol_order is a positive integer smaller than or equal to 16.
    if isfield(profile_options, ProfileOptionKey.MAX_TOL_ORDER.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)) || profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) <= 0 || profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) > 16
            error("MATLAB:checkValidityProfileOptions:max_tol_orderNotValid", "The field 'max_tol_order' of options should be a positive integer smaller than or equal to 16.");
        end
    end


    % Judge whether profile_options.max_eval_factor is a positive integer.
    if isfield(profile_options, ProfileOptionKey.MAX_EVAL_FACTOR.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value)) || profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) <= 0
            error("MATLAB:checkValidityProfileOptions:max_eval_factorNotValid", "The field 'max_eval_factor' of options should be a positive integer.");
        end
    end


    % Judge whether profile_options.merit_fun is a function handle.
    if isfield(profile_options, ProfileOptionKey.MERIT_FUN.value)
        if ~isa(profile_options.(ProfileOptionKey.MERIT_FUN.value), 'function_handle')
            error("MATLAB:checkValidityProfileOptions:merit_funNotValid", "The field 'merit_fun' of options should be a function handle.");
        end
    end


    % Judge whether profile_options.project_x0 is a boolean.
    if isfield(profile_options, ProfileOptionKey.PROJECT_X0.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.PROJECT_X0.value))
            error("MATLAB:checkValidityProfileOptions:project_x0NotValid", "The filed 'project_x0' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.run_plain is a boolean.
    if isfield(profile_options, ProfileOptionKey.RUN_PLAIN.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.RUN_PLAIN.value))
            error("MATLAB:checkValidityProfileOptions:run_plainNotValid", "The field 'run_plain' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.score_only is a boolean.
    if isfield(profile_options, ProfileOptionKey.SCORE_ONLY.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SCORE_ONLY.value))
            error("MATLAB:checkValidityProfileOptions:score_onlyNotValid", "The field 'score_only' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.summarize_performance_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_performance_profilesNotValid", "The field 'summarize_performance_profiles' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.summarize_data_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_data_profilesNotValid", "The field of 'summarize_data_profiles' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.summarize_log_ratio_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_log_ratio_profilesNotValid", "The field 'summarize_log_ratio_profiles' of options should be a boolean.");
        end
        if profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) && numel(solvers) > 2
            warning("MATLAB:checkValidityProfileOptions:summarize_log_ratio_profilesOnlyWhenTwoSolvers", "The log-ratio profiles are available only when there are exactly two solvers. We will not generate the log-ratio profiles.");
            profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) = false;
        end
    end


    % Judge whether profile_options.summarize_output_based_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_output_based_profilesNotValid", "The field 'summarize_output_based_profiles' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.silent is a boolean.
    if isfield(profile_options, ProfileOptionKey.SILENT.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SILENT.value))
            error("MATLAB:checkValidityProfileOptions:silentNotValid", "The field 'silent' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.solver_verbose is among 0, 1, and 2. If silent is true, solver_verbose should be 0 or 1 (only print errors). If it is 2, print a message saying that solver_verbose will be set to 1.
    if isfield(profile_options, ProfileOptionKey.SOLVER_VERBOSE.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value)) || ~ismember(profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value), [0, 1, 2])
            error("MATLAB:checkValidityProfileOptions:solver_verboseNotValid", "The field 'solver_verbose' of options should be either 0, 1, or 2.");
        end
        if profile_options.(ProfileOptionKey.SILENT.value) && profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
            warning("MATLAB:checkValidityProfileOptions:solver_verboseSetToZero", "solver_verbose will be set to 1 because silent is true.");
            profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) = 1;
        end
    end


    % Judge whether profile_options.semilogx is a boolean.
    if isfield(profile_options, ProfileOptionKey.SEMILOGX.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SEMILOGX.value))
            error("MATLAB:checkValidityProfileOptions:semilogxNotValid", "The field 'semilogx' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.normalized_scores is a boolean.
    if isfield(profile_options, ProfileOptionKey.NORMALIZED_SCORES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.NORMALIZED_SCORES.value))
            error("MATLAB:checkValidityProfileOptions:normalized_scoresNotValid", "The field 'normalized_scores' of options should be a boolean.");
        end
    end


    % Judge whether profile_options.score_weight_fun is a function handle.
    if isfield(profile_options, ProfileOptionKey.SCORE_WEIGHT_FUN.value)
        if ~isa(profile_options.(ProfileOptionKey.SCORE_WEIGHT_FUN.value), 'function_handle')
            error("MATLAB:checkValidityProfileOptions:score_weight_funNotValid", "The field 'score_weight_fun' of options should be a function handle.");
        end
    end


    % Judge whether profile_options.score_fun is a function handle.
    if isfield(profile_options, ProfileOptionKey.SCORE_FUN.value)
        if ~isa(profile_options.(ProfileOptionKey.SCORE_FUN.value), 'function_handle')
            error("MATLAB:checkValidityProfileOptions:score_funNotValid", "The field 'score_fun' of options should be a function handle.");
        end
    end


    % Judge whether profile_options.load is a char or a string.
    if isfield(profile_options, ProfileOptionKey.LOAD.value)
        if isempty(profile_options.(ProfileOptionKey.LOAD.value))
            profile_options.(ProfileOptionKey.LOAD.value) = '';
        end
        if ~ischarstr(profile_options.(ProfileOptionKey.LOAD.value))
            error("MATLAB:checkValidityProfileOptions:loadNotValid", "The field 'load' of options should be a char or a string.");
        end
        profile_options.(ProfileOptionKey.LOAD.value) = char(profile_options.(ProfileOptionKey.LOAD.value));
    end


    % Judge whether profile_options.solvers_to_load is a vector of different integers (with the same size as profile_options.solver_names if it exists).
    if isfield(profile_options, ProfileOptionKey.SOLVERS_TO_LOAD.value)
        if ~isfield(profile_options, ProfileOptionKey.LOAD.value) || isempty(profile_options.(ProfileOptionKey.LOAD.value))
            error("MATLAB:checkValidityProfileOptions:solvers_to_loadNoLoad", "The field 'solvers_to_load' of options can only be used when the field 'load' is set.");
        end
        profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value) = unique(profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value));
        if ~isintegervector(profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value)) || any(profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value) < 1) || numel(profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value)) < 2
            error("MATLAB:checkValidityProfileOptions:solvers_to_loadNotValid", "The field 'solvers_to_load' of options should be a vector of different integers greater than or equal to 1. At least two indices should be provided.");
        end
        if isfield(profile_options, ProfileOptionKey.SOLVER_NAMES.value)
            if numel(profile_options.(ProfileOptionKey.SOLVERS_TO_LOAD.value)) ~= numel(profile_options.(ProfileOptionKey.SOLVER_NAMES.value))
                error("MATLAB:checkValidityProfileOptions:solvers_to_loadAndsolver_namesLengthNotSame", "The field 'solvers_to_load' of options should have the same length as 'solver_names'.");
            end
        end
    end


    % Judge whether profile_options.line_colors is a cell array of strings or a matrix whose rows are RGB triplets.
    if isfield(profile_options, ProfileOptionKey.LINE_COLORS.value)
        % If it is a cell array of strings, check whether each string is a valid color short name.
        if iscell(profile_options.(ProfileOptionKey.LINE_COLORS.value))
            % Convert it to a row cell array of strings.
            if size(profile_options.(ProfileOptionKey.LINE_COLORS.value), 1) > 1
                profile_options.(ProfileOptionKey.LINE_COLORS.value) = profile_options.(ProfileOptionKey.LINE_COLORS.value)';
            end
            if ~all(cellfun(@(x) ischarstr(x) && ismember(x, {'r', 'g', 'b', 'c', 'm', 'y', 'k'}), profile_options.(ProfileOptionKey.LINE_COLORS.value)))
                error("MATLAB:checkValidityProfileOptions:line_colorsCellNotValid", "The field 'line_colors' of options can only contain 'r', 'g', 'b', 'c', 'm', 'y', 'k' when it is a cell array of strings.");
            end
            % Convert the cell array of strings to a matrix whose rows are RGB triplets.
            profile_options.(ProfileOptionKey.LINE_COLORS.value) = cellfun(@char, profile_options.(ProfileOptionKey.LINE_COLORS.value), 'UniformOutput', false);
            profile_options.(ProfileOptionKey.LINE_COLORS.value) = cellfun(@shortname2rgb, profile_options.(ProfileOptionKey.LINE_COLORS.value), 'UniformOutput', false);
            profile_options.(ProfileOptionKey.LINE_COLORS.value) = cell2mat(profile_options.(ProfileOptionKey.LINE_COLORS.value)');
        % If it is a matrix, check whether each row is an RGB triplet.
        elseif isnumeric(profile_options.(ProfileOptionKey.LINE_COLORS.value))
            if size(profile_options.(ProfileOptionKey.LINE_COLORS.value), 2) ~= 3 || ~all(all(profile_options.(ProfileOptionKey.LINE_COLORS.value) >= 0 & profile_options.(ProfileOptionKey.LINE_COLORS.value) <= 1))
                error("MATLAB:checkValidityProfileOptions:line_colorsMatrixNotValid", "The field 'line_colors' of options should be a matrix whose rows are RGB triplets when it is a matrix.");
            end
        else
            error("MATLAB:checkValidityProfileOptions:line_colorsNotValid", "The field 'line_colors' of options should be a cell array of short color names or a matrix whose rows are RGB triplets.");
        end
    end


    % Judge whether profile_options.line_styles is a cell array of chars or strings that are the combinations of line styles ('-', '-.', ':', '--') and markers ('o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'), where line styles cannot be 'none'.
    if isfield(profile_options, ProfileOptionKey.LINE_STYLES.value)
        % List all the valid combinations of line styles and markers.
        valid_line_style_marker_combinations = {...
        '-', '-o', '-+', '-*', '-.', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p', '-h', 'o-', '+-', '*-', '.-', 'x-', 's-', 'd-', '^-', 'v-', '>-', '<-', 'p-', 'h-',...
        '-.', '-.o', '-.+', '-.*', '-..', '-.x', '-.s', '-.d', '-.^', '-.v', '-.>', '-.<', '-.p', '-.h', 'o-.', '+-.', '*-.', '.-.', 'x-.', 's-.', 'd-.', '^-.', 'v-.', '>-.', '<-.', 'p-.', 'h-.',...
        ':', ':o', ':+', ':*', ':.', ':x', ':s', ':d', ':^', ':v', ':>', ':<', ':p', ':h', 'o:', '+:', '*:', '.:', 'x:', 's:', 'd:', '^:', 'v:', '>:', '<:', 'p:', 'h:',...
        '--', '--o', '--+', '--*', '--.', '--x', '--s', '--d', '--^', '--v', '-->', '--<', '--p', '--h', 'o--', '+--', '*--', '.--', 'x--', 's--', 'd--', '^--', 'v--', '>--', '<--', 'p--', 'h--'...
        };
        if ~iscell(profile_options.(ProfileOptionKey.LINE_STYLES.value)) || ~all(cellfun(@(x) ischarstr(x) && ismember(x, valid_line_style_marker_combinations), profile_options.(ProfileOptionKey.LINE_STYLES.value)))
            error("MATLAB:checkValidityProfileOptions:line_stylesNotValid", "The field 'line_styles' of options should be a cell array of chars or strings that are the combinations of line styles ('-', '-.', ':', '--') and markers ('none', 'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h').");
        end
        % Convert it to a row cell array of strings.
        if size(profile_options.(ProfileOptionKey.LINE_STYLES.value), 1) > 1
            profile_options.(ProfileOptionKey.LINE_STYLES.value) = profile_options.(ProfileOptionKey.LINE_STYLES.value)';
        end
    end


    % Judge whether profile_options.line_widths is a positive scalar or vector.
    if isfield(profile_options, ProfileOptionKey.LINE_WIDTHS.value)
        if ~isrealvector(profile_options.(ProfileOptionKey.LINE_WIDTHS.value)) || any(profile_options.(ProfileOptionKey.LINE_WIDTHS.value) <= 0)
            error("MATLAB:checkValidityProfileOptions:line_widthsNotValid", "The field 'line_widths' of options should be a positive scalar or vector.");
        end
    end


    % Judge whether profile_options.bar_colors is a cell array of strings or a matrix whose rows are RGB triplets.
    if isfield(profile_options, ProfileOptionKey.BAR_COLORS.value)
        % If it is a cell array of strings, check whether each string is a valid color short name.
        if iscell(profile_options.(ProfileOptionKey.BAR_COLORS.value))
            % Convert it to a row cell array of strings.
            if size(profile_options.(ProfileOptionKey.BAR_COLORS.value), 1) > 1
                profile_options.(ProfileOptionKey.BAR_COLORS.value) = profile_options.(ProfileOptionKey.BAR_COLORS.value)';
            end
            if ~all(cellfun(@(x) ischarstr(x) && ismember(x, {'r', 'g', 'b', 'c', 'm', 'y', 'k'}), profile_options.(ProfileOptionKey.BAR_COLORS.value)))
                error("MATLAB:checkValidityProfileOptions:bar_colorsCellNotValid", "The field 'bar_colors' of options can only contain 'r', 'g', 'b', 'c', 'm', 'y', 'k' when it is a cell array of strings.");
            end
            % Convert the cell array of strings to a matrix whose rows are RGB triplets.
            profile_options.(ProfileOptionKey.BAR_COLORS.value) = cellfun(@char, profile_options.(ProfileOptionKey.BAR_COLORS.value), 'UniformOutput', false);
            profile_options.(ProfileOptionKey.BAR_COLORS.value) = cellfun(@shortname2rgb, profile_options.(ProfileOptionKey.BAR_COLORS.value), 'UniformOutput', false);
            profile_options.(ProfileOptionKey.BAR_COLORS.value) = cell2mat(profile_options.(ProfileOptionKey.BAR_COLORS.value)');
        % If it is a matrix, check whether each row is an RGB triplet.
        elseif isnumeric(profile_options.(ProfileOptionKey.BAR_COLORS.value))
            if size(profile_options.(ProfileOptionKey.BAR_COLORS.value), 2) ~= 3 || ~all(all(profile_options.(ProfileOptionKey.BAR_COLORS.value) >= 0 & profile_options.(ProfileOptionKey.BAR_COLORS.value) <= 1))
                error("MATLAB:checkValidityProfileOptions:bar_colorsMatrixNotValid", "The field 'bar_colors' of options should be a matrix whose rows are RGB triplets when it is a matrix.");
            end
        else
            error("MATLAB:checkValidityProfileOptions:bar_colorsNotValid", "The field 'bar_colors' of options should be a cell array of short color names or a matrix whose rows are RGB triplets.");
        end
    end
end

function rgb = shortname2rgb(shortname)
%SHORTNAME2RGB Convert a short color name to an RGB triplet

    switch shortname
        case 'r'
            rgb = [1, 0, 0];
        case 'g'
            rgb = [0, 1, 0];
        case 'b'
            rgb = [0, 0, 1];
        case 'c'
            rgb = [0, 1, 1];
        case 'm'
            rgb = [1, 0, 1];
        case 'y'
            rgb = [1, 1, 0];
        case 'k'
            rgb = [0, 0, 0];
    end
end