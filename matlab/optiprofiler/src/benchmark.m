function benchmark(solvers, varargin)
%BENCHMARK Create multiple profiles for benchmarking optimization solvers on a
%   set of problems with different features.
%
%   BENCHMARK(SOLVERS) creates performance profiles and data profiles for the
%   given SOLVERS with default unconstrained problem set and feature plain.
%
%   BENCHMARK(SOLVERS, FEATURE_NAMES) creates performance profiles and data
%   profiles for the given SOLVERS with default unconstrained problem set and
%   specified features FEATURE_NAMES.
%
%   BENCHMARK(SOLVERS, FEATURE_NAMES, PROBLEM) creates history plots for the
%   given SOLVERS with the specified one PROBLEM and features FEATURE_NAMES.
%
%   BENCHMARK(SOLVERS, OPTIONS) creates performance profiles and data profiles
%   for the given SOLVERS with options specified in the struct OPTIONS.
%   OPTIONS can contain the following fields:
%       1. options for profiles and plots:
%       - n_jobs: the number of parallel jobs to run the test. Default is 1.
%       - benchmark_id: the identifier of the test. It is used to create the
%         specific directory to store the results. Default is '.' .
%       - range_type: the type of the uncertainty interval. For stochastic
%         features, we run several times of the experiments and get average
%         curves and uncertainty intervals. Default is 'minmax', meaning that we
%         takes the pointwise minimum and maximum of the curves.
%       - std_factor: the factor multiplied to the standard deviation to get the
%         uncertainty interval in the case where range_type is 'std', meaning
%         that we take the pointwise mean plus/minus the standard deviation of
%         the curves. Default is 1.
%       - savepath: the path to store the results. Default is the current
%         directory where the function is called.
%       - max_tol_order: the maximum order of the tolerance. In any profile (in
%         our case, performance profiles, data profiles, and log-ratio profiles),
%         we need to set a group of 'tolerances' to define the 'convergence' of
%         the solvers. (Details can be found in the references.) We will set the
%         tolerances as 10^(-1:-1:-max_tol_order). Default is 10.
%       - max_eval_factor: the factor multiplied to each problem's dimension to
%         get the maximum number of evaluations for each problem. Default is 500.
%       - project_x0: whether to project the initial point to the feasible set.
%         Default is false.
%       - run_plain: whether to run the plain feature when the feature is
%         stochastic (e.g., the feature is 'noisy' and you set the 'run_plain'
%         to true, then we will additionally run the 'plain' feature and use the
%         results to define the 'convergence' of the solvers). Default is true.
%       - summarize_performance_profiles: whether to add all the performance
%         profiles to the summary PDF. Default is true.
%       - summarize_data_profiles: whether to add all the data profiles to the
%         summary PDF. Default is true.
%       - summarize_log_ratio_profiles: whether to add all the log-ratio profiles
%         to the summary PDF. Default is false.
%       - summarize_output_based_profiles: whether to add all the output-based
%         profiles of the selected profiles to the summary PDF. Default is true.
%       - summarize_funhist: whether to add the history plots of the objective
%         function values to the summary PDF. Default is true. (Note: it is only
%         available when there is only one problem.)
%       - summarize_maxcvhist: whether to add the history plots of the maximum
%         constraint violation to the summary PDF. Default is true. (Note: it is
%         only available when there is only one problem.)
%       - summarize_merithist: whether to add the history plots of the merit
%         values to the summary PDF. Default is true. (Note: it is only available
%         when there is only one problem.)
%       - summarize_cumminhist: whether to add all the cumulative minimum of the 
%         selected history plots to the summary PDF. Default is true. (Note: it
%         is only available when there is only one problem.)
%       2. options for features:
%       - n_runs: the number of runs of the experiments under the given feature.
%         Default is 10 for stochastic features and 1 for deterministic features.
%       - distribution: the distribution of random vectors in stochastic
%         features. It should be a function handle,
%               (random stream, dimension) -> random vector)
%         accepting a random stream and the dimension of a problem, and returning
%         a vector with the same dimension. Default is the standard multivariate
%         normal distribution.
%       - noise_level: the magnitude of the noise in stochastic features. Default
%         is 10^-3.
%       - noise_type: the type of the noise in stochastic features. It should be
%         either 'absolute' or 'relative'. Default is 'relative'.
%       - significant_digits: the number of significant digits in the 'truncated'
%         feature. Default is 6.
%       - perturbed_trailing_zeros: whether we will set the trailing zeros of the
%         objective function value to be perturbed (randomly generated) in the
%         'perturbed_x0' feature. Default is true.
%       - rotated: whether to use a random or given rotation matrix to rotate the
%         coordinates of a problem in the 'linearly_transformed' feature. Default
%         is true.
%       - invertible_transformation: the invertible transformation in the
%         'linearly_transformed' feature. It should be a function handle,
%               (random stream, dimension) -> (invertible matrix, inverse matrix)
%         accepting a random stream and the dimension of a problem, and returning
%         an invertible matrix and its inverse. Default is generating a random
%         orthogonal matrix following the uniform distribution on SO(dimension),
%         and the inverse is the transpose of the matrix.
%       - condition_number: the condition number of a scaling matrix, which will
%         be composed with the objective function in the 'linearly_transformed'
%         feature. It should be a function handle,
%               dimension -> condition_number
%         accepting the dimension of a problem and returning a scalar. Default is
%         @(n) 1, meaning that the scaling matrix is the identity matrix.
%       - unrelaxable_bounds: whether the bound constraints are unrelaxable or
%         not in the 'unrelaxable_constraints' feature. Default is false.
%       - unrelaxable_linear_constraints: whether the linear constraints are
%         unrelaxable or not in the 'unrelaxable_constraints' feature. Default is
%         false.
%       - unrelaxable_nonlinear_constraints: whether the nonlinear constraints
%         are unrelaxable or not in the 'unrelaxable_constraints' feature.
%         Default is false.
%       - rate_nan: the probability that the evaluation of the objective function
%         will return NaN in the 'random_nan' feature. Default is 0.05.
%       - modifier: the modifier function to modify the objective function value
%         in the 'custom' feature. It should be a function handle,
%               (current point, current objective function value, seed) ->
%               modified objective function value
%         accepting the current point, the current objective function value, and
%         a random seed, and returning the modified objective function value. No
%         default setting.
%       3. options for CUTEst:
%       Note that the CUTEst we used is the MATLAB codes from a GitHub repository
%       called 'S2MPJ', created by Professor Serge Gratton and Professor Philippe
%       L. Toint. More details can be found in the following website.
%           https://github.com/GrattonToint/S2MPJ
%       - problem_type: the type of the problems to be selected. It should be a
%         string containing the combination of 'u' (unconstrained), 'b' (bound
%         constrained), 'l' (linearly constrained), and 'n' (nonlinearly
%         constrained). Default is 'u'.
%       - mindim: the minimum dimension of the problems to be selected. Default
%         is 1.
%       - maxdim: the maximum dimension of the problems to be selected. Default
%         is 5.
%       - mincon: the minimum number of constraints of the problems to be
%         selected. Default is 0.
%       - maxcon: the maximum number of constraints of the problems to be
%         selected. Default is 10 for linearly or nonlinearly constrained
%         problems, and 0 for the others.
%       - excludelist: the list of problems to be excluded. Default is not to
%         exclude any problem.
%
%   For more information of performance and data profiles, see [1]_, [2]_, [4]_. 
%   For that of log-ratio profiles, see [3]_, [5]_.
%   Pay attention that log-ratio profiles are available only when there are
%   exactly two solvers.
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
%
%   **************************************************************************
%   Authors:    Cunxin HUANG (cun-xin.huang@connect.polyu.hk)
%               Tom M. RAGONNEAU (t.ragonneau@gmail.com)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%   **************************************************************************



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the input arguments. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin == 0
        error("MATLAB:benchmark:solverMustBeProvided", "solvers must be provided.");
    elseif nargin == 1
        % When input only contains one argument, we assume the user chooses benchmark(solvers) and test plain feature.
        feature_names = 'plain';
        labels = cellfun(@func2str, solvers, 'UniformOutput', false);
        cutest_problem_names = {};
        custom_problem_loader = {};
        custom_problem_names = {};
        options = struct();
    elseif nargin == 2
        if ischarstr(varargin{1}) || (iscell(varargin{1}) && all(cellfun(@ischarstr, varargin{1})))
            % When input contains two arguments and the second argument is a char or cell of char, we assume the user chooses benchmark(solvers, feature_names).
            feature_names = varargin{1};
            labels = cellfun(@func2str, solvers, 'UniformOutput', false);
            cutest_problem_names = {};
            custom_problem_loader = {};
            custom_problem_names = {};
            options = struct();
        elseif isstruct(varargin{1})
            % When input contains two arguments and the second argument is a struct, we assume the user chooses benchmark(solvers, options).
            options = varargin{1};
            if isfield(options, 'feature_names')
                feature_names = options.feature_names;
                options = rmfield(options, 'feature_names');
            else
                feature_names = 'plain';
            end
            if isfield(options, 'labels')
                labels = options.labels;
                options = rmfield(options, 'labels');
            else
                labels = cellfun(@func2str, solvers, 'UniformOutput', false);
            end
            if isfield(options, 'cutest_problem_names')
                cutest_problem_names = options.cutest_problem_names;
                options = rmfield(options, 'cutest_problem_names');
            else
                cutest_problem_names = {};
            end
            if isfield(options, 'custom_problem_loader') && isfield(options, 'custom_problem_names')
                custom_problem_loader = options.custom_problem_loader;
                custom_problem_names = options.custom_problem_names;
                options = rmfield(options, 'custom_problem_loader');
                options = rmfield(options, 'custom_problem_names');
            elseif isfield(options, 'custom_problem_loader') || isfield(options, 'custom_problem_names')
                error("MATLAB:benchmark:LoaderAndNamesNotSameTime", "custom_problem_loader and custom_problem_names must be provided at the same time.");
            else
                custom_problem_loader = {};
                custom_problem_names = {};
            end
        else
            error("MATLAB:benchmark:SecondArgumentWrongType", "The second argument must be a cell array of feature names or a struct of options.");
        end
    elseif nargin == 3
        % When input contains three arguments, we assume the user chooses benchmark(solvers, feature_names, problem).
        if ~isa(varargin{2}, 'Problem')
            error("MATLAB:benchmark:ThirdArgumentNotProblem", "The third argument must be a Problem object.");
        end
        feature_names = varargin{1};
        labels = cellfun(@func2str, solvers, 'UniformOutput', false);
        cutest_problem_names = {};
        custom_problem_loader = @(x) varargin{2};
        custom_problem_names = {'custom'};
        options = struct('n_jobs', 1);
    else
        error("MATLAB:benchmark:TooMuchInput", "Invalid number of arguments. The function must be called with one, two, or three arguments.");
    end

    % Preprocess the solvers.
    if ~iscell(solvers) || ~all(cellfun(@(s) isa(s, 'function_handle'), solvers))
        error("MATLAB:benchmark:solversWrongType", "The solvers must be a cell array of function handles.");
    end
    if numel(solvers) < 2
        error("MATLAB:benchmark:solversAtLeastTwo", "At least two solvers must be given.");
    end

    % Preprocess the feature_names.
    if ~ischarstr(feature_names) && ~(iscell(feature_names) && all(cellfun(@ischarstr, feature_names)))
        % feature_names must be a char or string, or a cell array of chars or strings.
        error("MATLAB:benchmark:feature_namesNotcharstrOrCellOfcharstr", "The feature names must be a cell array of chars or strings.");
    end
    if ischarstr(feature_names)
        % Convert the char or string to a cell array of chars.
        feature_names = {char(feature_names)};
    end
    if strcmp(feature_names, 'all')
        % If 'all' is given, add all the feature names except 'custom'.
        feature_names = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
        custom_idx = strcmp(feature_names, FeatureName.CUSTOM.value);
        feature_names(custom_idx) = [];
    end
    if ~all(cellfun(@(x) ismember(x, {enumeration('FeatureName').value}), feature_names))
        error("MATLAB:benchmark:feature_namesNotValid", "The feature names must be valid.");
    end
    % Remove the duplicates.
    feature_names = unique(feature_names, 'stable');

    % Preprocess the labels.
    if ~iscell(labels) || ~all(cellfun(@(l) ischarstr(l), labels))
        error("MATLAB:benchmark:labelsNotCellOfcharstr", "The labels must be a cell of chars or strings.");
    end
    if numel(labels) ~= 0 && numel(labels) ~= numel(solvers)
        error("MATLAB:benchmark:labelsAndsolversLengthNotSame", "The number of labels must equal the number of solvers.");
    end
    if numel(labels) == 0
        labels = cellfun(@func2str, solvers, 'UniformOutput', false);
    end

    % Preprocess the custom problems.
    if ~isempty(custom_problem_loader) && ~isa(custom_problem_loader, 'function_handle')
        error("MATLAB:benchmark:customloaderNotFunctionHandle", "The custom problem loader must be a function handle.");
    end
    if ~isempty(custom_problem_loader)
        if isempty(custom_problem_names)
            error("MATLAB:benchmark:customnamesCanNotBeEmptyWhenHavingcustomloader", "The custom problem names must be provided.");
        else
            try
                [~, p] = evalc('custom_problem_loader(custom_problem_names{1})');
            catch
                p = [];
            end
            if isempty(p) || ~isa(p, 'Problem')
                error("MATLAB:benchmark:customloaderNotAcceptcustomnames", "The custom problem loader must be able to accept one signature 'custom_problem_names'. The first problem %s could not be loaded, or custom problem loader did not return a Problem object.", custom_problem_names{1});
            end
        end
    elseif ~isempty(custom_problem_names)
        error("MATLAB:benchmark:customloaderCanNotBeEmptyWhenHavingcustomnames", "A custom problem loader must be given to load custom problems.");
    end
    if ~isempty(custom_problem_names)
        if ~ischarstr(custom_problem_names) && ~(iscell(custom_problem_names) && all(cellfun(@ischarstr, custom_problem_names)))
            error("MATLAB:benchmark:customnamesNotcharstrOrCellOfcharstr", "The custom problem names must be a cell array of chars or strings.");
        end
        if ischarstr(custom_problem_names)
            custom_problem_names = {custom_problem_names};
        end
        custom_problem_names = cellfun(@char, custom_problem_names, 'UniformOutput', false);  % Convert to cell array of chars.
    end

    % Set default profile options.
    profile_options = getDefaultProfileOptions();
    
    % Initialize the options for feature and cutest.
    feature_options = struct();
    cutest_options = struct();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Parse options for feature, cutest, and profile. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fieldNames = fieldnames(options);
    for i_field = 1:numel(fieldNames)
        key = fieldNames{i_field};
        value = options.(key);

        % We can also use: validFeatureOptionKeys = {enumeration('FeatureOptionKey').value};  
        % However, it only works for MATLAB R2021b or later.
        validFeatureOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('FeatureOptionKey')), 'UniformOutput', false);
        validCutestOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('CutestOptionKey')), 'UniformOutput', false);
        validProfileOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('ProfileOptionKey')), 'UniformOutput', false);

        if ismember(key, validFeatureOptionKeys)
            feature_options.(key) = value;
        elseif ismember(key, validCutestOptionKeys)
            cutest_options.(key) = value;
        elseif ismember(key, validProfileOptionKeys)
            profile_options.(key) = value;
        else
            error("MATLAB:benchmark:UnknownOptions", "Unknown option: %s", key);
        end
    end

    % Note: the validity of the feature options has been checked in the Feature constructor, so we do not need to check it here.

    % Check whether the cutest options are valid.
    checkValidityCutestOptions(cutest_options);

    % Check whether the profile options are valid.
    checkValidityProfileOptions(profile_options, solvers);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Use cutest_options to select problems. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cutest_problem_names = loadCutestNames(cutest_options, cutest_problem_names, custom_problem_loader, custom_problem_names);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set the default values for plotting. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create the directory to store the results.
    path_out = setSavingPath(profile_options);

    % Set the default values for plotting.
    set(groot, 'DefaultLineLineWidth', 1);
    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start the computation of the profiles. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note: if user gives more than one feature, we only allow to benchmark under all the default feature options. (We do not want to mix the feature options.)
    if length(feature_names) > 1 && numel(fieldnames(feature_options)) > 0
        error("MATLAB:benchmark:OnlyOneFeatureWhenHavingfeature_options", "Only one feature can be specified when feature options are given.");
    end

    for i_feature = 1:length(feature_names)
        feature_name = feature_names{i_feature};

        % Build feature.
        feature = Feature(feature_name, feature_options);
        fprintf('INFO: Starting the computation of the "%s" profiles.\n', feature.name);

        % Paths to store the results.
        path_feature = fullfile(path_out, feature.name);
        if ~exist(path_feature, 'dir')
            mkdir(path_feature);
        end

        % Solve all the problems.
        [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, time_processes] = solveAllProblems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, profile_options);
        merit_histories = computeMeritValues(fun_histories, maxcv_histories, maxcv_init);
        merit_out = computeMeritValues(fun_out, maxcv_out, maxcv_init);
        merit_init = computeMeritValues(fun_init, maxcv_init, maxcv_init);

        % If there are no problems solved, skip the rest of the code, print a message, and continue with the next feature.
        if isempty(problem_names)
            fprintf('INFO: No problems were solved for the "%s" feature.\n', feature.name);
            continue;
        end

        % If there is only one problem given, draw profiles of fun_histories, maxcv_histories, and merit_histories.
        if length(cutest_problem_names) + length(custom_problem_names) == 1
            % Sqeeze the dimension of the 'problem' axis.
            fun_histories = reshape(fun_histories, size(fun_histories, 2), size(fun_histories, 3), size(fun_histories, 4));
            maxcv_histories = reshape(maxcv_histories, size(maxcv_histories, 2), size(maxcv_histories, 3), size(maxcv_histories, 4));
            merit_histories = reshape(merit_histories, size(merit_histories, 2), size(merit_histories, 3), size(merit_histories, 4));

            % Create the figure for the summary.
            warning('off');
            n_cols = 0;
            is_fun = profile_options.(ProfileOptionKey.SUMMARIZE_FUNHIST.value);
            is_maxcv = profile_options.(ProfileOptionKey.SUMMARIZE_MAXCVHIST.value);
            is_merit = profile_options.(ProfileOptionKey.SUMMARIZE_MERITHIST.value);
            is_cum = profile_options.(ProfileOptionKey.SUMMARIZE_CUMMINHIST.value);
            if is_fun
                n_cols = n_cols + 1;
            end
            if is_maxcv
                n_cols = n_cols + 1;
            end
            if is_merit
                n_cols = n_cols + 1;
            end
            if is_cum
                multiplier = 2;
            else
                multiplier = 1;
            end
            defaultFigurePosition = get(0, 'DefaultFigurePosition');
            default_width = defaultFigurePosition(3);
            default_height = defaultFigurePosition(4);
            fig_summary = figure('Position', [defaultFigurePosition(1:2), n_cols * default_width, multiplier * default_height], 'visible', 'off');
            T_summary = tiledlayout(fig_summary, multiplier, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
            T_title = strrep(feature.name, '_', '\_');
            P_title = strrep(problem_names{1}, '_', '\_');
            title(T_summary, ['Solving ``', P_title, '" with the ``', T_title, '" feature'], 'Interpreter', 'latex', 'FontSize', 14);
            % Use gobjects to create arrays of handles and axes.
            t_summary = gobjects(multiplier, 1);
            axs_summary = gobjects([multiplier, 1, 1, n_cols]);
            i_axs = 0;
            for i = 1:multiplier
                t_summary(i) = tiledlayout(T_summary, 1, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
                t_summary(i).Layout.Tile = i;
                for j = 1:n_cols
                    i_axs = i_axs + 1;
                    axs_summary(i_axs) = nexttile(t_summary(i));
                end
            end
            ylabel(t_summary(1), "History profiles", 'Interpreter', 'latex', 'FontSize', 14);
            if is_cum
                ylabel(t_summary(2), "Cummin history profiles", 'Interpreter', 'latex', 'FontSize', 14);
            end

            cell_axs_summary_cum = {};
            if is_fun && is_maxcv && is_merit
                cell_axs_summary = {axs_summary(1), axs_summary(2), axs_summary(3)};
                if is_cum
                    cell_axs_summary_cum = {axs_summary(4), axs_summary(5), axs_summary(6)};
                end
            elseif (is_fun && is_maxcv) || (is_fun && is_merit) || (is_maxcv && is_merit)
                cell_axs_summary = {axs_summary(1), axs_summary(2)};
                if is_cum
                    cell_axs_summary_cum = {axs_summary(3), axs_summary(4)};
                end
            elseif is_fun || is_maxcv || is_merit
                cell_axs_summary = {axs_summary(1)};
                if is_cum
                    cell_axs_summary_cum = {axs_summary(2)};
                end
            end

            pdf_summary = fullfile(path_out, 'summary.pdf');

            [fig_fun, fig_maxcv, fig_merit] = drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, labels, cell_axs_summary, true, is_fun, is_maxcv, is_merit, false, profile_options);
            [fig_cummin_fun, fig_cummin_maxcv, fig_cummin_merit] = drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, labels, cell_axs_summary_cum, is_cum, is_fun, is_maxcv, is_merit, true, profile_options);

            eps_fun = fullfile(path_feature, 'fun_hist.eps');
            print(fig_fun, eps_fun, '-depsc');
            pdf_fun = fullfile(path_feature, 'fun_hist.pdf');
            print(fig_fun, pdf_fun, '-dpdf');
            eps_maxcv = fullfile(path_feature, 'maxcv_hist.eps');
            print(fig_maxcv, eps_maxcv, '-depsc');
            pdf_maxcv = fullfile(path_feature, 'maxcv_hist.pdf');
            print(fig_maxcv, pdf_maxcv, '-dpdf');
            eps_merit = fullfile(path_feature, 'merit_hist.eps');
            print(fig_merit, eps_merit, '-depsc');
            pdf_merit = fullfile(path_feature, 'merit_hist.pdf');
            print(fig_merit, pdf_merit, '-dpdf');
            eps_fun_cum = fullfile(path_feature, 'cummin_fun_hist.eps');
            print(fig_cummin_fun, eps_fun_cum, '-depsc');
            pdf_fun_cum = fullfile(path_feature, 'cummin_fun_hist.pdf');
            print(fig_cummin_fun, pdf_fun_cum, '-dpdf');
            eps_maxcv_cum = fullfile(path_feature, 'cummin_maxcv_hist.eps');
            print(fig_cummin_maxcv, eps_maxcv_cum, '-depsc');
            pdf_maxcv_cum = fullfile(path_feature, 'cummin_maxcv_hist.pdf');
            print(fig_cummin_maxcv, pdf_maxcv_cum, '-dpdf');
            eps_merit_cum = fullfile(path_feature, 'cummin_merit_hist.eps');
            print(fig_cummin_merit, eps_merit_cum, '-depsc');
            pdf_merit_cum = fullfile(path_feature, 'cummin_merit_hist.pdf');
            print(fig_cummin_merit, pdf_merit_cum, '-dpdf');

            close(fig_fun);
            close(fig_maxcv);
            close(fig_merit);
            close(fig_cummin_fun);
            close(fig_cummin_maxcv);
            close(fig_cummin_merit);

            if i_feature == 1
                exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector');
            else
                exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector', 'Append', true);
            end
            fprintf('Detailed results stored in %s\n', path_feature);
    
            warning('on');
            continue;

        end

        % Determine the least merit value for each problem.
        merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
        if feature.isStochastic && profile_options.(ProfileOptionKey.RUN_PLAIN.value)
            feature_plain = Feature(FeatureName.PLAIN.value);
            fprintf('INFO: Starting the computation of the "plain" profiles.\n');
            [fun_histories_plain, maxcv_histories_plain, ~, ~, ~, ~, ~, ~, ~, time_processes_plain] = solveAllProblems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature_plain, profile_options);
            time_processes = time_processes + time_processes_plain;
            merit_histories_plain = computeMeritValues(fun_histories_plain, maxcv_histories_plain, maxcv_init);
            merit_min_plain = min(min(min(merit_histories_plain, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
            merit_min = min(merit_min, merit_min_plain, 'omitnan');
        end

        % Create the directories for the performance profiles, data profiles, and log-ratio profiles.
        path_perf = fullfile(path_feature, 'figures', 'perf');
        path_data = fullfile(path_feature, 'figures', 'data');
        path_log_ratio = fullfile(path_feature, 'figures', 'log-ratio');
        path_perf_hist = fullfile(path_perf, 'history-based');
        path_data_hist = fullfile(path_data, 'history-based');
        path_log_ratio_hist = fullfile(path_log_ratio, 'history-based');
        path_perf_out = fullfile(path_perf, 'output-based');
        path_data_out = fullfile(path_data, 'output-based');
        path_log_ratio_out = fullfile(path_log_ratio, 'output-based');
        if ~exist(path_perf_hist, 'dir')
            mkdir(path_perf_hist);
        end
        if ~exist(path_data_hist, 'dir')
            mkdir(path_data_hist);
        end
        if ~exist(path_perf_out, 'dir')
            mkdir(path_perf_out);
        end
        if ~exist(path_data_out, 'dir')
            mkdir(path_data_out);
        end

        % Store the names of the problems.
        path_txt = fullfile(path_feature, 'problems.txt');
        [~, idx] = sort(lower(problem_names));
        sorted_problem_names = problem_names(idx);
        sorted_time_processes = time_processes(idx);
        fid = fopen(path_txt, 'w');
        if fid == -1
            error("MATLAB:benchmark:FileCannotOpen", "Cannot open the file %s.", path_txt);
        end
        for i = 1:length(sorted_problem_names)
            count = fprintf(fid, "%s: %.2f seconds\n", sorted_problem_names{i}, sorted_time_processes(i));
            if count < 0
                error("MATLAB:benchmark:FailToEditFile", "Failed to record data for %s.", sorted_problem_names{i});
            end
        end
        fclose(fid);

        [n_problems, n_solvers, n_runs, ~] = size(merit_histories);
        if n_solvers <= 2
            if ~exist(path_log_ratio_hist, 'dir')
                mkdir(path_log_ratio_hist);
            end
            if ~exist(path_log_ratio_out, 'dir')
                mkdir(path_log_ratio_out);
            end
        end

        max_tol_order = profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value);
        tolerances = 10.^(-1:-1:-max_tol_order);
        pdf_summary = fullfile(path_out, 'summary.pdf');
        pdf_perf_hist_summary = fullfile(path_perf, 'perf_hist.pdf');
        pdf_perf_out_summary = fullfile(path_perf, 'perf_out.pdf');
        pdf_data_hist_summary = fullfile(path_data, 'data_hist.pdf');
        pdf_data_out_summary = fullfile(path_data, 'data_out.pdf');
        pdf_log_ratio_hist_summary = fullfile(path_log_ratio, 'log-ratio_hist.pdf');
        pdf_log_ratio_out_summary = fullfile(path_log_ratio, 'log-ratio_out.pdf');

        % Create the figure for the summary.
        warning('off');
        n_rows = 0;
        is_perf = profile_options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value);
        is_data = profile_options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value);
        is_log_ratio = profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) && (n_solvers <= 2);
        is_output_based = profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value);
        if is_perf
            n_rows = n_rows + 1;
        end
        if is_data
            n_rows = n_rows + 1;
        end
        if is_log_ratio
            n_rows = n_rows + 1;
        end
        if is_output_based
            multiplier = 2;
        else
            multiplier = 1;
        end
        defaultFigurePosition = get(0, 'DefaultFigurePosition');
        default_width = defaultFigurePosition(3);
        default_height = defaultFigurePosition(4);
        fig_summary = figure('Position', [defaultFigurePosition(1:2), profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * default_width, multiplier * n_rows * default_height], 'visible', 'off');
        T_summary = tiledlayout(fig_summary, multiplier, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        T_title = strrep(feature.name, '_', '\_');
        title(T_summary, ['Profiles with the ``', T_title, '" feature'], 'Interpreter', 'latex', 'FontSize', 24);
        % Use gobjects to create arrays of handles and axes.
        t_summary = gobjects(multiplier, 1);
        axs_summary = gobjects([multiplier, 1, n_rows, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)]);
        i_axs = 0;
        for i = 1:multiplier
            t_summary(i) = tiledlayout(T_summary, n_rows, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value), 'Padding', 'compact', 'TileSpacing', 'compact');
            t_summary(i).Layout.Tile = i;
            for j = 1:n_rows * profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)
                i_axs = i_axs + 1;
                axs_summary(i_axs) = nexttile(t_summary(i));
            end
        end
        ylabel(t_summary(1), "History-based profiles", 'Interpreter', 'latex', 'FontSize', 20);
        if is_output_based
            ylabel(t_summary(2), "Output-based profiles", 'Interpreter', 'latex', 'FontSize', 20);
        end

        for i_profile = 1:profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)
            tolerance = tolerances(i_profile);
            [tolerance_str, tolerance_latex] = formatFloatScientificLatex(tolerance);
            fprintf("Creating profiles for tolerance %s.\n", tolerance_str);
            tolerance_label = ['$\mathrm{tol} = ' tolerance_latex '$'];

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
            cell_axs_summary_out = {};
            if is_perf && is_data && is_log_ratio
                cell_axs_summary_hist = {axs_summary(i_profile), axs_summary(i_profile + max_tol_order), axs_summary(i_profile + 2 * max_tol_order)};
                if is_output_based
                    cell_axs_summary_out = {axs_summary(i_profile + 3 * max_tol_order), axs_summary(i_profile + 4 * max_tol_order), axs_summary(i_profile + 5 * max_tol_order)};
                end
            elseif (is_perf && is_data) || (is_perf && is_log_ratio) || (is_data && is_log_ratio)
                cell_axs_summary_hist = {axs_summary(i_profile), axs_summary(i_profile + max_tol_order)};
                if is_output_based
                    cell_axs_summary_out = {axs_summary(i_profile + 2 * max_tol_order), axs_summary(i_profile + 3 * max_tol_order)};
                end
            elseif is_perf || is_data || is_log_ratio
                cell_axs_summary_hist = {axs_summary(i_profile)};
                if is_output_based
                    cell_axs_summary_out = {axs_summary(i_profile + max_tol_order)};
                end
            end

            [fig_perf_hist, fig_data_hist, fig_log_ratio_hist] = drawProfiles(work_hist, problem_dimensions, labels, tolerance_label, cell_axs_summary_hist, true, is_perf, is_data, is_log_ratio, profile_options);
            [fig_perf_out, fig_data_out, fig_log_ratio_out] = drawProfiles(work_out, problem_dimensions, labels, tolerance_label, cell_axs_summary_out, is_output_based, is_perf, is_data, is_log_ratio, profile_options);
            eps_perf_hist = fullfile(path_perf_hist, ['perf_hist_' int2str(i_profile) '.eps']);
            print(fig_perf_hist, eps_perf_hist, '-depsc');
            pdf_perf_hist = fullfile(path_perf_hist, ['perf_hist_' int2str(i_profile) '.pdf']);
            print(fig_perf_hist, pdf_perf_hist, '-dpdf');
            eps_perf_out = fullfile(path_perf_out, ['perf_out_' int2str(i_profile) '.eps']);
            print(fig_perf_out, eps_perf_out, '-depsc');
            pdf_perf_out = fullfile(path_perf_out, ['perf_out_' int2str(i_profile) '.pdf']);
            print(fig_perf_out, pdf_perf_out, '-dpdf');
            eps_data_hist = fullfile(path_data_hist, ['data_hist_' int2str(i_profile) '.eps']);
            print(fig_data_hist, eps_data_hist, '-depsc');
            pdf_data_hist = fullfile(path_data_hist, ['data_hist_' int2str(i_profile) '.pdf']);
            print(fig_data_hist, pdf_data_hist, '-dpdf');
            eps_data_out = fullfile(path_data_out, ['data_out_' int2str(i_profile) '.eps']);
            print(fig_data_out, eps_data_out, '-depsc');
            pdf_data_out = fullfile(path_data_out, ['data_out_' int2str(i_profile) '.pdf']);
            print(fig_data_out, pdf_data_out, '-dpdf');
            if n_solvers <= 2
                eps_log_ratio_hist = fullfile(path_log_ratio_hist, ['log-ratio_hist_' int2str(i_profile) '.eps']);
                print(fig_log_ratio_hist, eps_log_ratio_hist, '-depsc');
                pdf_log_ratio_hist = fullfile(path_log_ratio_hist, ['log-ratio_hist_' int2str(i_profile) '.pdf']);
                print(fig_log_ratio_hist, pdf_log_ratio_hist, '-dpdf');
            end
            if n_solvers <= 2
                eps_log_ratio_out = fullfile(path_log_ratio_out, ['log-ratio_out_' int2str(i_profile) '.eps']);
                print(fig_log_ratio_out, eps_log_ratio_out, '-depsc');
                pdf_log_ratio_out = fullfile(path_log_ratio_out, ['log-ratio_out_' int2str(i_profile) '.pdf']);
                print(fig_log_ratio_out, pdf_log_ratio_out, '-dpdf');
            end
            if i_profile == 1
                exportgraphics(fig_perf_hist, pdf_perf_hist_summary, 'ContentType', 'vector');
                exportgraphics(fig_perf_out, pdf_perf_out_summary, 'ContentType', 'vector');
                exportgraphics(fig_data_hist, pdf_data_hist_summary, 'ContentType', 'vector');
                exportgraphics(fig_data_out, pdf_data_out_summary, 'ContentType', 'vector');
                if n_solvers <= 2
                    exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist_summary, 'ContentType', 'vector');
                end
                if ~isempty(fig_log_ratio_out)
                    exportgraphics(fig_log_ratio_out, pdf_log_ratio_out_summary, 'ContentType', 'vector');
                end
            else
                exportgraphics(fig_perf_hist, pdf_perf_hist_summary, 'ContentType', 'vector', 'Append', true);
                exportgraphics(fig_perf_out, pdf_perf_out_summary, 'ContentType', 'vector', 'Append', true);
                exportgraphics(fig_data_hist, pdf_data_hist_summary, 'ContentType', 'vector', 'Append', true);
                exportgraphics(fig_data_out, pdf_data_out_summary, 'ContentType', 'vector', 'Append', true);
                if n_solvers <= 2
                    exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist_summary, 'ContentType', 'vector', 'Append', true);
                end
                if ~isempty(fig_log_ratio_out)
                    exportgraphics(fig_log_ratio_out, pdf_log_ratio_out_summary, 'ContentType', 'vector', 'Append', true);
                end
            end

            % Close the figures.
            close(fig_perf_hist);
            close(fig_perf_out);
            close(fig_data_hist);
            close(fig_data_out);
            if n_solvers <= 2
                close(fig_log_ratio_hist);
            end
            if ~isempty(fig_log_ratio_out)
                close(fig_log_ratio_out);
            end

        end

        if i_feature == 1
            exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector');
        else
            exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector', 'Append', true);
        end
        fprintf('Detailed results stored in %s\n', path_feature);

        warning('on');

    end

    % Close the figures.
    close(fig_summary);
    fprintf('Summary stored in %s\n', path_out);

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