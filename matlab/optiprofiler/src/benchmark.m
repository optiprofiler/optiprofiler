function [solver_scores, profile_scores, problem_scores, profiles] = benchmark(solvers, varargin)
%BENCHMARK Create multiple profiles for benchmarking optimization solvers on a
%   set of problems with different features.
%
%   Signatures:
%
%   SOLVER_SCORES = BENCHMARK(SOLVERS) creates performance profiles, data
%   profiles, and log-ratio profiles (later say profiles) for the given SOLVERS
%   on the default unconstrained problem set, returning SOLVER_SCORES based on
%   the profiles. SOLVERS is a cell array of function handles. We require
%   SOLVERS to accept specified inputs and return specified outputs. Details
%   can be found in the following 'Cautions' part.
%
%   SOLVER_SCORES = BENCHMARK(SOLVERS, FEATURE_NAME) creates profiles and data
%   profiles for the given SOLVERS on the default unconstrained problem set
%   with the specified feature FEATURE_NAME.
%
%   SOLVER_SCORES = BENCHMARK(SOLVERS, OPTIONS) creates profiles for the given
%   SOLVERS with options specified in the struct OPTIONS. See 'Options' part
%   for more details.
%
%   [SOLVER_SCORES, PROFILE_SCORES] = BENCHMARK(...) returns a 4D tensor
%   PROFILE_SCORES containing scores for all profiles. See `scoring_fun` in
%   'Options' part for more details.
%
%   [SOLVER_SCORES, PROFILE_SCORES, PROBLEM_SCORES] = BENCHMARK(...) returns a
%   3D tensor PROBLEM_SCORES containing scores of the solvers on all problems.
%
%   [SOLVER_SCORES, PROFILE_SCORES, PROBLEM_SCORES, PROFILES] = BENCHMARK(...)
%   returns a cell array PROFILES containing all the figure objects of the
%   profiles.
%
%   Options:
%
%   Options should be specified in a struct. The following are the available
%   fields of the struct:
%
%       1. options for profiles and plots:
%       - n_jobs: the number of parallel jobs to run the test. Default is 1.
%       - benchmark_id: the identifier of the test. It is used to create the
%         specific directory to store the results. Default is '.' .
%       - feature_stamp: the stamp of the feature with the given options. It is
%         used to create the specific directory to store the results. Default
%         is different for different features.
%       - range_type: the type of the uncertainty interval. For stochastic
%         features, we run several times of the experiments and get average
%         curves and uncertainty intervals. Default is 'minmax', meaning that
%         we takes the pointwise minimum and maximum of the curves.
%       - savepath: the path to store the results. Default is the current
%         directory where the function is called.
%       - max_tol_order: the maximum order of the tolerance. In any profile
%         (performance profiles, data profiles, and log-ratio profiles), we
%         need to set a group of 'tolerances' to define the 'convergence' of
%         the solvers. (Details can be found in the references.) We will set
%         the tolerances as `10^(-1:-1:-max_tol_order)`. Default is 10.
%       - max_eval_factor: the factor multiplied to each problem's dimension to
%         get the maximum number of evaluations for each problem. Default is
%         500.
%       - project_x0: whether to project the initial point to the feasible set.
%         Default is false.
%       - run_plain: whether to run an extra experiment with the 'plain'
%         feature. Default is false.
%       - draw_plots: whether to draw history plots and profiles. Default is
%         true.
%       - summarize_performance_profiles: whether to add all the performance
%         profiles to the summary PDF. Default is true.
%       - summarize_data_profiles: whether to add all the data profiles to the
%         summary PDF. Default is true.
%       - summarize_log_ratio_profiles: whether to add all the log-ratio
%         profiles to the summary PDF. Default is false.
%       - summarize_output_based_profiles: whether to add all the output-based
%         profiles of the selected profiles to the summary PDF. Default is
%         true.
%       - semilogx: whether to use the semilogx scale during calculating the
%         integral of the performance profiles and data profiles. Default is
%         true.
%       - scoring_fun: the scoring function to calculate the scores of the
%         solvers. It should be a function handle as follows:
%               profile_scores -> solver_scores,
%         where `profile_scores` is a 4D tensor containing scores for all
%         profiles. The first dimension is the index of the solver, the second
%         is the index of tolerance starting from 1, the third represents
%         history-based or output-based profiles, and the fourth represents
%         performance profiles, data profiles, or log-ratio profiles. The
%         default scoring function takes the average of the history-based
%         performance profiles along the tolerance axis and then normalizes the
%         average by the maximum value of the average.
%
%       2. options for features:
%       - feature_name: the name of the feature. Default is 'plain'.
%       - n_runs: the number of runs of the experiments under the given
%         feature. Default is 10 for stochastic features and 1 for
%         deterministic features.
%       - distribution: the distribution of perturbation in 'perturbed_x0'
%         feature or noise in 'noisy' feature. It should be either a string
%         (or char), or a function handle
%               (random stream, dimension) -> random vector,
%         accepting a random stream and the dimension of a problem and
%         returning a random vector with the given dimension. In 'perturbed_x0'
%         case, the char should be either 'spherical' or 'gaussian' (default is
%         'spherical'). In 'noisy' case, the char should be either 'gaussian'
%         or 'uniform' (default is 'gaussian').
%       - noise_level: the magnitude of the noise in stochastic features.
%         Default is 10^-3.
%       - noise_type: the type of the noise in stochastic features. It should
%         be either 'absolute', 'relative', or 'mixed'. Default is 'mixed'.
%       - significant_digits: the number of significant digits in the
%         'truncated' feature. Default is 6.
%       - perturbed_trailing_zeros: whether we will randomize the trailing
%         zeros of the objective function value in the 'perturbed_x0' feature.
%         Default is false.
%       - rotated: whether to use a random or given rotation matrix to rotate
%         the coordinates of a problem in the 'linearly_transformed' feature.
%         Default is true.
%       - condition_factor: the scaling factor of the condition number of the
%         linear transformation in the 'linearly_transformed' feature. More
%         specifically, the condition number of the linear transformation will
%         2 ^ (condition_factor * n / 2), where `n` is the dimension of the
%         problem. Default is 0.
%       - nan_rate: the probability that the evaluation of the objective
%         function will return NaN in the 'random_nan' feature. Default is
%         0.05.
%       - unrelaxable_bounds: whether the bound constraints are unrelaxable or
%         not in the 'unrelaxable_constraints' feature. Default is false.
%       - unrelaxable_linear_constraints: whether the linear constraints are
%         unrelaxable or not in the 'unrelaxable_constraints' feature. Default
%         is false.
%       - unrelaxable_nonlinear_constraints: whether the nonlinear constraints
%         are unrelaxable or not in the 'unrelaxable_constraints' feature.
%         Default is false.
%       - mesh_size: the size of the mesh in the 'quantized' feature. Default
%         is 10^-3.
%       - ground_truth: whether the feature is the ground truth or not. Default
%         is true.
%       - mod_x0: the modifier function to modify the inital guess in the 
%         'custom' feature. It should be a function handle as follows:
%               (random_stream, problem) -> modified_x0,
%         where `problem` is an instance of the class Problem, and
%         `modified_x0` is the modified initial guess. No default.
%       - mod_affine: the modifier function to generate the affine
%         transformation applied to the variables in the 'custom' feature. It
%         should be a function handle as follows:
%               (random_stream, problem) -> (A, b, inv),
%         where `problem` is an instance of the class Problem, `A` is the
%         matrix of the affine transformation, `b` is the vector of the affine
%         transformation, and `inv` is the inverse of matrix `A`. No default.
%       - mod_bounds: the modifier function to modify the bound constraints in
%         the 'custom' feature. It should be a function handle as follows:
%               (random_stream, problem) -> (modified_xl, modified_xu),
%         where `problem` is an instance of the class Problem, `modified_xl` is
%         the modified lower bound, and `modified_xu` is the modified upper
%         bound. No default.
%       - mod_linear_ub: the modifier function to modify the linear inequality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               (random_stream, problem) -> (modified_aub, modified_bub),
%         where `problem` is an instance of the class Problem, `modified_aub`
%         is the modified matrix of the linear inequality constraints, and
%         `modified_bub` is the modified vector of the linear inequality
%         constraints. No default.
%       - mod_linear_eq: the modifier function to modify the linear equality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               (random_stream, problem) -> (modified_aeq, modified_beq),
%         where `problem` is an instance of the class Problem, `modified_aeq`
%         is the modified matrix of the linear equality constraints, and
%         `modified_beq` is the modified vector of the linear equality
%         constraints. No default.
%       - mod_fun: the modifier function to modify the objective function in
%         the 'custom' feature. It should be a function handle as follows:
%               (x, random_stream, problem) -> modified_fun,
%         where `x` is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_fun` is the modified objective function
%         value. No default.
%       - mod_cub: the modifier function to modify the nonlinear inequality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               (x, random_stream, problem) -> modified_cub,
%         where x is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_cub` is the modified vector of the
%         nonlinear inequality constraints. No default.
%       - mod_ceq: the modifier function to modify the nonlinear equality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               (x, random_stream, problem) -> modified_ceq,
%         where x is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_ceq` is the modified vector of the
%         nonlinear equality constraints. No default.
%
%       3. options for CUTEst:
%       Note that the CUTEst we used is the MATLAB codes from a GitHub
%       repository called 'S2MPJ', created by Professor Serge Gratton and
%       Professor Philippe L. Toint. More details can be found in the following
%       website.
%           https://github.com/GrattonToint/S2MPJ
%       - p_type: the type of the problems to be selected. It should be a
%         string containing the combination of 'u' (unconstrained), 'b' (bound
%         constrained), 'l' (linearly constrained), and 'n' (nonlinearly
%         constrained). Default is 'u'.
%       - mindim: the minimum dimension of the problems to be selected. Default
%         is 1.
%       - maxdim: the maximum dimension of the problems to be selected. Default
%         is 2.
%       - mincon: the minimum number of constraints of the problems to be
%         selected. Default is 0.
%       - maxcon: the maximum number of constraints of the problems to be
%         selected. Default is 10 for linearly or nonlinearly constrained
%         problems, and 0 for the others.
%       - excludelist: the list of problems to be excluded. Default is not to
%         exclude any problem.
%
%       4. other options:
%       - solver_names: the names of the solvers. Default is the function names
%         of the solvers.
%       - solver_isrand: whether the solvers are randomized or not. Default is
%         a logical array of the same length as the number of solvers, where
%         the value is true if the solver is randomized, and false otherwise.
%       - problem: a instance of the class Problem. If it is provided, we will
%         only solve this problem and generate the history plots for it.
%         Default is not to provide any problem.
%       - cutest_problem_names: the names of the problems in the CUTEst library
%         to be selected. Default is not to select any problem from the CUTEst
%         library by name but by the options above.
%       - custom_problem_loader: the function handle to load the custom
%         problems. It should be a function handle as follows:
%               (problem_name) -> problem,
%         where `problem_name` is the name of the problem, and `problem` is an
%         instance of the class Problem. Default is not to load any custom
%         problem.
%       - custom_problem_names: the names of the custom problems to be
%         selected. Default is not to select any custom problem.
%
%   Cautions:
%
%   1. Each solver in SOLVERS should accept the following signature(s):
%       - for an unconstrained problem,
%           x = solver(fun, x0),
%         where `fun` is a function handle of the objective function accepting
%         a column vector and returning a real number, and `x0` is the initial
%         guess which is a column vector;
%       - for a bound-constrained problem,
%           x = solver(fun, x0, xl, xu),
%         where `xl` and `xu` are the lower and upper bounds of the variables
%         which are column vectors (they can contain Inf or -Inf);
%       - for a linearly constrained problem,
%           x = solver(fun, x0, xl, xu, aub, bub, aeq, beq);
%         where `aub` and `aeq` are the matrices of the linear inequality and
%         equality constraints, and `bub` and `beq` are the vectors of the
%         linear inequality and equality constraints;
%       - for a nonlinearly constrained problem,
%           x = solver(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq),
%         where `cub` and `ceq` are the functions of the nonlinear inequality
%         and equality constraints accepting a column vector and returning a
%         column vector.
%   2. The log-ratio profiles are available only when there are exactly two 
%      solvers.
%
%   For more information of performance and data profiles, see [1]_, [2]_,
%   [4]_. For that of log-ratio profiles, see [3]_, [5]_.
%
%   References:
%   .. [1] E. D. Dolan and J. J. Moré. Benchmarking optimization software with
%          performance profiles. *Math. Program.*, 91(2):201–213, 2002.
%          `doi:10.1007/s101070100263
%          <https://doi.org/10.1007/s101070100263>.
%   .. [2] N. Gould and J. Scott. A note on performance profiles for
%          benchmarking software. *ACM Trans. Math. Software*, 43(2):15:1–5,
%          2016. `doi:10.1145/2950048 <https://doi.org/10.1145/2950048>.
%   .. [3] S. Gratton and Ph. L. Toint. S2MPJ and CUTEst optimization problems
%          for Matlab, Python and Julia. arXiv:2407.07812, 2024.
%   .. [4] J. L. Morales. A numerical study of limited memory BFGS methods.
%          *Appl. Math. Lett.*, 15(4):481–487, 2002.
%          `doi:10.1016/S0893-9659(01)00162-8
%          <https://doi.org/10.1016/S0893-9659(01)00162-8>.
%   .. [5] J. J. Moré and S. M. Wild. Benchmarking derivative-free optimization
%          algorithms. *SIAM J. Optim.*, 20(1):172–191, 2009.
%          `doi:10.1137/080724083 <https://doi.org/10.1137/080724083>.
%   .. [6] H.-J. M. Shi, M. Q. Xuan, F. Oztoprak, and J. Nocedal. On the
%          numerical performance of finite-difference-based methods for
%          derivative-free optimization. *Optim. Methods Softw.*,
%          38(2):289–311, 2023. `doi:10.1080/10556788.2022.2121832
%          <https://doi.org/10.1080/10556788.2022.2121832>.
%
%   **************************************************************************
%   Authors:    Cunxin HUANG (cun-xin.huang@connect.polyu.hk)
%               Tom M. RAGONNEAU (t.ragonneau@gmail.com)
%               Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%   **************************************************************************

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preprocess the input arguments. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 0
        error("MATLAB:benchmark:solverMustBeProvided", "A cell of function handles (callable solvers) must be provided as the first argument.");
    elseif nargin == 1
        % When input only contains one argument, we assume the user chooses benchmark(solvers) and
        % test plain feature.
        feature_name = 'plain';
        options = struct();
    elseif nargin == 2
        if ischarstr(varargin{1}) || (iscell(varargin{1}) && all(cellfun(@ischarstr, varargin{1})))
            % When input contains two arguments and the second argument is a char or cell of char,
            % we assume the user chooses benchmark(solvers, feature_name).
            feature_name = varargin{1};
            options = struct();
        elseif isstruct(varargin{1})
            % When input contains two arguments and the second argument is a struct, we assume the
            % user chooses benchmark(solvers, options).
            options = varargin{1};
            options_store = options;
            if isfield(options, 'feature_name')
                feature_name = options.feature_name;
                options = rmfield(options, 'feature_name');
            else
                feature_name = 'plain';
            end
        else
            error("MATLAB:benchmark:SecondArgumentWrongType", ...
            "The second argument for `benchmark` must be a feature name or a struct of options.");
        end
    else
        error("MATLAB:benchmark:TooMuchInput", ...
        "Invalid number of arguments. `benchmark` function at most takes two arguments.");
    end

    % Preprocess the solvers.
    if ~iscell(solvers) || ~all(cellfun(@(s) isa(s, 'function_handle'), solvers))
        error("MATLAB:benchmark:solversWrongType", "The first argument for `benchmark` must be a cell array of function handles.");
    end
    if numel(solvers) < 2
        error("MATLAB:benchmark:solversAtLeastTwo", "The first argument for `benchmark` must be a cell array of at least two function handles since we need to compare at least two solvers.");
    end

    % Preprocess the feature_name.
    if ~ischarstr(feature_name)
        % feature_name must be a char or string.
        error("MATLAB:benchmark:feature_nameNotcharstr", "The field `feature_name` of `options` for `benchmark` must be a char or string.");
    end
    feature_name = char(lower(feature_name));
    valid_feature_names = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
    if ~ismember(feature_name, valid_feature_names)
        error("MATLAB:benchmark:feature_nameNotValid", "The field `feature_name` of `options` for `benchmark` must be one of the valid feature names: %s.", strjoin(valid_feature_names, ', '));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Parse options for feature, cutest, profile and others %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    feature_options = struct();
    cutest_options = struct();
    profile_options = struct();
    other_options = struct();

    fieldNames = fieldnames(options);
    for i_field = 1:numel(fieldNames)
        key = fieldNames{i_field};
        value = options.(key);

        % We can also use: validFeatureOptionKeys = {enumeration('FeatureOptionKey').value};  
        % However, it only works for MATLAB R2021b or later.
        validFeatureOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('FeatureOptionKey')), 'UniformOutput', false);
        validCutestOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('CutestOptionKey')), 'UniformOutput', false);
        validProfileOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('ProfileOptionKey')), 'UniformOutput', false);
        validOtherOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('OtherOptionKey')), 'UniformOutput', false);

        if ismember(key, validFeatureOptionKeys)
            feature_options.(key) = value;
        elseif ismember(key, validCutestOptionKeys)
            cutest_options.(key) = value;
        elseif ismember(key, validProfileOptionKeys)
            profile_options.(key) = value;
        elseif ismember(key, validOtherOptionKeys)
            other_options.(key) = value;
        else
            error("MATLAB:benchmark:UnknownOptions", "Unknown `option` for `benchmark`: %s", key);
        end
    end

    % Build feature.
    feature = Feature(feature_name, feature_options);

    % Set default values for the unspecified options.
    other_options = getDefaultOtherOptions(solvers, other_options);
    cutest_options = getDefaultCutestOptions(cutest_options, other_options);
    profile_options = getDefaultProfileOptions(feature, profile_options);

    % Check the validity of the options.
    % Note: the validity of the feature options has been checked in the Feature constructor, so we
    % do not need to check it here.
    cutest_options = checkValidityCutestOptions(cutest_options);
    profile_options = checkValidityProfileOptions(solvers, profile_options);
    other_options = checkValidityOtherOptions(solvers, other_options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize output variables. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_solvers = numel(solvers);
    solver_scores = zeros(n_solvers, 1);
    profile_scores = [];
    problem_scores = [];
    profiles = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Use cutest_options to select problems. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = loadCutestNames(cutest_options, other_options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Set the default values for plotting. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create the directory to store the results.
    path_out = setSavingPath(profile_options);

    % Set the default values for plotting.
    set(groot, 'DefaultLineLineWidth', 1);
    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Create the directory to store the results. %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create the directory to store the results. If it already exists, overwrite it.
    path_feature = fullfile(path_out, profile_options.(ProfileOptionKey.FEATURE_STAMP.value));
    if ~exist(path_feature, 'dir')
        mkdir(path_feature);
    else
        rmdir(path_feature, 's');
        mkdir(path_feature);
    end
    if ~profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
        path_hist_plots = '';
    else
        path_hist_plots = fullfile(path_feature, 'history_plots');
        if ~exist(path_hist_plots, 'dir')
            mkdir(path_hist_plots);
        end
    end

    % Create the directory to store options and log files.
    path_log = fullfile(path_feature, 'test_log');
    if ~exist(path_log, 'dir')
        mkdir(path_log);
    else
        rmdir(path_log, 's');
        mkdir(path_log);
    end
    try
        if exist('options_store', 'var')
            save(fullfile(path_log, 'options_store.mat'), 'options_store');
        end
    catch
        fprintf("INFO: Failed to save the `options` of the current experiment.\n");
    end
    log_file = fullfile(path_log, 'log.txt');
    diary(log_file);

    % If `n_runs` is not specified, the feature is deterministic, and at least one solver is
    % randomized, then we set `n_runs` to 5.
    if ~isfield(feature_options, FeatureOptionKey.N_RUNS.value) && ~feature.is_stochastic && any(other_options.(OtherOptionKey.SOLVER_ISRAND.value))
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf("INFO: We set `n_runs` to 5 since the feature is deterministic and at least one solver is randomized and `n_runs` is not specified.\n\n");
        end
        feature.options.(FeatureOptionKey.N_RUNS.value) = 5;
    end

    % We try to copy the script or function that calls the benchmark function to the log directory.
    try
        calling_script = dbstack(1, '-completenames');
        if ~isempty(calling_script)
            copyfile(calling_script.file, path_log);
            fprintf("INFO: The script or function that calls `benchmark` function is copied to: %s.\n\n", path_log);
        end
    catch
        fprintf("INFO: Failed to copy the script or function that calls `benchmark` function to the log directory.\n\n");
    end

    if profile_options.(ProfileOptionKey.DRAW_PLOTS.value) && ~isfield(other_options, OtherOptionKey.PROBLEM.value)
        % Create the directories for the performance profiles, data profiles, and log-ratio profiles.
        path_perf_hist = fullfile(path_feature, 'detailed_profiles', 'perf_history-based');
        path_data_hist = fullfile(path_feature, 'detailed_profiles', 'data_history-based');
        path_log_ratio_hist = fullfile(path_feature, 'detailed_profiles', 'log-ratio_history-based');
        path_perf_out = fullfile(path_feature, 'detailed_profiles', 'perf_output-based');
        path_data_out = fullfile(path_feature, 'detailed_profiles', 'data_output-based');
        path_log_ratio_out = fullfile(path_feature, 'detailed_profiles', 'log-ratio_output-based');
        
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
        if length(solvers) == 2
            if ~exist(path_log_ratio_hist, 'dir')
                mkdir(path_log_ratio_hist);
            end
            if ~exist(path_log_ratio_out, 'dir')
                mkdir(path_log_ratio_out);
            end
        end

        pdf_perf_hist_summary = fullfile(path_feature, 'perf_hist.pdf');
        pdf_perf_out_summary = fullfile(path_feature, 'perf_out.pdf');
        pdf_data_hist_summary = fullfile(path_feature, 'data_hist.pdf');
        pdf_data_out_summary = fullfile(path_feature, 'data_out.pdf');
        pdf_log_ratio_hist_summary = fullfile(path_feature, 'log-ratio_hist.pdf');
        pdf_log_ratio_out_summary = fullfile(path_feature, 'log-ratio_out.pdf');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solve all the problems. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('INFO: Starting the computation of the "%s" profiles.\n', feature.name);
    end
    [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, time_processes, problem_unsolved] = solveAllProblems(solvers, feature, profile_options, other_options, true, path_hist_plots);
    merit_histories = computeMeritValues(fun_histories, maxcv_histories, maxcv_init);
    merit_out = computeMeritValues(fun_out, maxcv_out, maxcv_init);
    merit_init = computeMeritValues(fun_init, maxcv_init, maxcv_init);

    % If there are no problems solved, skip the rest of the code, print a message, and return.
    if isempty(problem_names)
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('INFO: No problems were solved for the "%s" feature.\n', feature.name);
        end
        diary off;
        return;
    end

    % If a specific problem is provided to `other_options`, we only solve this problem and generate
    % the history plots for it.
    if isfield(other_options, OtherOptionKey.PROBLEM.value)
        if profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
            % We move the history plots to the feature directory.
            try
                movefile(fullfile(path_hist_plots, '*'), path_feature);
                rmdir(path_hist_plots, 's');
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fprintf('\nINFO: Detailed results stored in %s\n', path_feature);
                end
            catch
            end
        end

        % Compute `solver_scores` based on the `merit_histories` and `merit_init`.
        % We first find the least merit value for each problem. Then we set the score of the solver
        % having the least merit value to 1, and the score of the solver having the largest merit
        % value to 0. The scores of the other solvers are linearly interpolated between 0 and 1.
        solver_merit_mins = squeeze(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'));
        solver_merit_mins_min = min(solver_merit_mins, [], 2, 'omitnan');
        solver_scores = (merit_init - solver_merit_mins) ./ (merit_init - solver_merit_mins_min);
        solver_scores = solver_scores';
        diary off;
        return;
    end

    % Determine the least merit value for each problem.
    merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
    if feature.run_plain && profile_options.(ProfileOptionKey.RUN_PLAIN.value)
        feature_plain = Feature(FeatureName.PLAIN.value);
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('\nINFO: Starting the computation of the profiles under "plain" feature.\n');
        end
        [fun_histories_plain, maxcv_histories_plain, ~, ~, ~, ~, ~, problem_names_plain, ~, time_processes_plain] = solveAllProblems(solvers, feature_plain, profile_options, other_options, false, {});
        merit_histories_plain = computeMeritValues(fun_histories_plain, maxcv_histories_plain, maxcv_init);
        merit_min_plain = min(min(min(merit_histories_plain, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
        for i_problem = 1:numel(problem_names)
            % Check whether the problem is solved under the plain feature.
            if ismember(problem_names{i_problem}, problem_names_plain)
                % Find the index of the problem in the plain feature.
                idx = find(strcmp(problem_names{i_problem}, problem_names_plain), 1);
                merit_min(i_problem) = min(merit_min(i_problem), merit_min_plain(idx), 'omitnan');
                time_processes(i_problem) = time_processes(i_problem) + time_processes_plain(idx);
            end
        end
    end

    % Store the names of the problems.
    path_txt = fullfile(path_log, 'problems.txt');
    [~, idx] = sort(lower(problem_names));
    sorted_problem_names = problem_names(idx);
    sorted_time_processes = time_processes(idx);
    max_name_length = max(cellfun(@length, sorted_problem_names));
    sorted_time_processes = num2cell(sorted_time_processes);
    max_time_length = max(cellfun(@(x) length(sprintf('%.2f', x)), sorted_time_processes));
    try
        fid = fopen(path_txt, 'w');
        for i = 1:length(sorted_problem_names)
            count = fprintf(fid, "%-*s      %*s\n", max_name_length, sorted_problem_names{i}, max_time_length, sprintf('%.2f seconds', sorted_time_processes{i}));
            if count < 0
                fprintf("INFO: Failed to record data for %s.", sorted_problem_names{i});
            end
        end
        if ~isempty(problem_unsolved)
            fprintf(fid, "\n");
            fprintf(fid, "Unsolved problems:\n");
            for i = 1:length(problem_unsolved)
                count = fprintf(fid, "%s\n", problem_unsolved{i});
                if count < 0
                    fprintf("INFO: Failed to record data for %s.", problem_unsolved{i});
                end
            end
        end
        fclose(fid);
    catch
        fprintf("INFO: Error occurred when writing the problem names to %s.\n", path_txt);
    end

    if profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
        % Merge all the pdf files in path_hist_plots to a single pdf file.
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('\nINFO: Merging all the history plots to a single PDF file.\n');
        end
        try
            mergePdfs(path_hist_plots, 'history_plots_summary.pdf', path_feature);
        catch
            warning('INFO: Failed to merge the history plots to a single PDF file.');
        end

        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('\nINFO: Detailed results stored in %s\n\n', path_feature);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Start the computation of all profiles. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [n_problems, n_solvers, n_runs, ~] = size(merit_histories);
    max_tol_order = profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value);
    tolerances = 10.^(-1:-1:-max_tol_order);

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

    if n_rows > 0
        fig_summary = figure('Position', [defaultFigurePosition(1:2), profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * default_width, multiplier * n_rows * default_height], 'visible', 'off');
        T_summary = tiledlayout(fig_summary, multiplier, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        T_feature_stamp = strrep(profile_options.(ProfileOptionKey.FEATURE_STAMP.value), '_', '\_');
        T_title = ['Profiles with the ``', T_feature_stamp, '" feature'];
        summary_width = profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * default_width;
        summary_height = multiplier * n_rows * default_height;
        summary_fontsize = min(summary_width / 75 * 2, summary_width / profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * 3 / 75 * 2);
        title(T_summary, T_title, 'Interpreter', 'latex', 'FontSize', summary_fontsize);
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
        ylabel(t_summary(1), "History-based profiles", 'Interpreter', 'latex', 'FontSize', min(summary_height / 60 * 2, summary_height / n_rows * 2 / 60 * 2));
        if is_output_based
            ylabel(t_summary(2), "Output-based profiles", 'Interpreter', 'latex', 'FontSize', min(summary_height / 60 * 2, summary_height / n_rows * 2 / 60 * 2));
        end
    end

    % Store the curves of the performance profiles, data profiles, and log-ratio profiles.
    profiles_perf = cell(2, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));
    profiles_data = cell(2, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));
    profiles_log_ratio = cell(2, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));

    curves = cell(1, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('\n');
    end

    for i_tol = 1:profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)
        curves{i_tol} = struct();
        curves{i_tol}.hist = struct();
        curves{i_tol}.hist.perf = cell(n_solvers, n_runs + 1);
        curves{i_tol}.hist.data = cell(n_solvers, n_runs + 1);
        if n_solvers == 2
            curves{i_tol}.hist.log_ratio = cell(1, 2);
        end
        curves{i_tol}.out = curves{i_tol}.hist;

        tolerance = tolerances(i_tol);
        [tolerance_str, tolerance_latex] = formatFloatScientificLatex(tolerance, 1);

        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf("INFO: Creating profiles for tolerance %s.\n", tolerance_str);
        end
        tolerance_label = ['$\mathrm{tol} = ' tolerance_latex '$'];

        work_hist = NaN(n_problems, n_solvers, n_runs);
        work_out = NaN(n_problems, n_solvers, n_runs);
        for i_problem = 1:n_problems
            for i_solver = 1:n_solvers
                for i_run = 1:n_runs
                    if isfinite(merit_min(i_problem))
                        threshold = max(tolerance * merit_init(i_problem) + (1 - tolerance) * merit_min(i_problem), ...
                        merit_min(i_problem));
                    else
                        threshold = -Inf;
                    end
                    if min(merit_histories(i_problem, i_solver, i_run, :), [], 'omitnan') <= threshold
                        work_hist(i_problem, i_solver, i_run) = find(merit_histories(i_problem, i_solver, i_run, :) ...
                        <= threshold, 1, 'first');
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
            cell_axs_summary_hist = {axs_summary(i_tol), axs_summary(i_tol + max_tol_order), ...
            axs_summary(i_tol + 2 * max_tol_order)};
            if is_output_based
                cell_axs_summary_out = {axs_summary(i_tol + 3 * max_tol_order), ...
                axs_summary(i_tol + 4 * max_tol_order), axs_summary(i_tol + 5 * max_tol_order)};
            end
        elseif (is_perf && is_data) || (is_perf && is_log_ratio) || (is_data && is_log_ratio)
            cell_axs_summary_hist = {axs_summary(i_tol), axs_summary(i_tol + max_tol_order)};
            if is_output_based
                cell_axs_summary_out = {axs_summary(i_tol + 2 * max_tol_order), ...
                axs_summary(i_tol + 3 * max_tol_order)};
            end
        elseif is_perf || is_data || is_log_ratio
            cell_axs_summary_hist = {axs_summary(i_tol)};
            if is_output_based
                cell_axs_summary_out = {axs_summary(i_tol + max_tol_order)};
            end
        else
            cell_axs_summary_hist = {};
            cell_axs_summary_out = {};
        end

        processed_solver_names = cellfun(@(s) strrep(s, '_', '\_'), other_options.(OtherOptionKey.SOLVER_NAMES.value), 'UniformOutput', false);
        [fig_perf_hist, fig_data_hist, fig_log_ratio_hist, curves{i_tol}.hist] = drawProfiles(work_hist, problem_dimensions, processed_solver_names, tolerance_label, cell_axs_summary_hist, true, is_perf, is_data, is_log_ratio, profile_options, curves{i_tol}.hist);
        [fig_perf_out, fig_data_out, fig_log_ratio_out, curves{i_tol}.out] = drawProfiles(work_out, problem_dimensions, processed_solver_names, tolerance_label, cell_axs_summary_out, is_output_based, is_perf, is_data, is_log_ratio, profile_options, curves{i_tol}.out);

        profiles_perf{1, i_tol} = fig_perf_hist;
        profiles_perf{2, i_tol} = fig_perf_out;
        profiles_data{1, i_tol} = fig_data_hist;
        profiles_data{2, i_tol} = fig_data_out;
        profiles_log_ratio{1, i_tol} = fig_log_ratio_hist;
        profiles_log_ratio{2, i_tol} = fig_log_ratio_out;

        if profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
            pdf_perf_hist = fullfile(path_perf_hist, ['perf_hist_', int2str(i_tol), '.pdf']);
            exportgraphics(fig_perf_hist, pdf_perf_hist, 'ContentType', 'vector');
            pdf_perf_out = fullfile(path_perf_out, ['perf_out_', int2str(i_tol), '.pdf']);
            exportgraphics(fig_perf_out, pdf_perf_out, 'ContentType', 'vector');
            pdf_data_hist = fullfile(path_data_hist, ['data_hist_', int2str(i_tol), '.pdf']);
            exportgraphics(fig_data_hist, pdf_data_hist, 'ContentType', 'vector');
            pdf_data_out = fullfile(path_data_out, ['data_out_', int2str(i_tol), '.pdf']);
            exportgraphics(fig_data_out, pdf_data_out, 'ContentType', 'vector');
            if n_solvers <= 2
                pdf_log_ratio_hist = fullfile(path_log_ratio_hist, ['log-ratio_hist_', int2str(i_tol), '.pdf']);
                exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist, 'ContentType', 'vector');
                pdf_log_ratio_out = fullfile(path_log_ratio_out, ['log-ratio_out_', int2str(i_tol), '.pdf']);
                exportgraphics(fig_log_ratio_out, pdf_log_ratio_out, 'ContentType', 'vector');
            end
            if i_tol == 1
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
        end

        % Close the figures.
        close(fig_perf_hist);
        close(fig_perf_out);
        close(fig_data_hist);
        close(fig_data_out);
        if n_solvers == 2
            close(fig_log_ratio_hist);
        end
        if ~isempty(fig_log_ratio_out)
            close(fig_log_ratio_out);
        end
    end

    profiles = cell(1, 4);
    profiles{1} = profiles_perf;
    profiles{2} = profiles_data;
    profiles{3} = profiles_log_ratio;
    if n_rows > 0
        % Check the operating system. If it is Windows, we will adjust the position, papersize, paperposition of the summary figure!
        if ispc
            fig_summary.Position = [defaultFigurePosition(1:2), profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * default_width, multiplier * n_rows * default_height];
            fig_summary.Units = 'centimeters';
            fig_summary.PaperUnits = 'centimeters';
            fig_summary.PaperSize = fig_summary.Position(3:4);
            fig_summary.PaperPosition = [0, 0, fig_summary.Position(3:4)];
            fig_summary.Children.Title.FontSize = min(fig_summary.Position(3) / 75 * 10, fig_summary.Position(3) / profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * 3 / 75 * 10);
            fig_summary.Children.Children(1).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
            fig_summary.Children.Children(2).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
        end
        profiles{4} = fig_summary;
    else
        profiles{4} = [];
    end

    % Store `profiles` in a mat file in the path_log directory.
    save(fullfile(path_log, 'profiles.mat'), 'profiles');

    % Compute the `solver_scores`.
    profile_scores = computeScores(curves, profile_options.(ProfileOptionKey.SEMILOGX.value));
    scoring_fun = profile_options.(ProfileOptionKey.SCORING_FUN.value);
    solver_scores = scoring_fun(profile_scores);
    save(fullfile(path_log, 'profile_scores.mat'), 'profile_scores');

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        % Print the scores of the solvers.
        fprintf('\n');
        fprintf('INFO: Scores of the solvers:\n');
        max_solver_name_length = max(cellfun(@length, processed_solver_names));
        for i_solver = 1:n_solvers
            format_info_str = sprintf('INFO: %%-%ds:    %%.4f\n', max_solver_name_length);
            fprintf(format_info_str, processed_solver_names{i_solver}, solver_scores(i_solver));
        end
    end

    if profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
        % Store the summary pdf. We will name the summary pdf as "summary_feature_name.pdf" and store it under
        % path_feature. We will also put a "summary.pdf" in the path_out directory, which will be a merged pdf of all
        % the "summary_feature_name.pdf" under path_out following the order of the feature_stamp.
        if n_rows > 0
            if ispc
                print(fig_summary, fullfile(path_feature, ['summary_' feature_name '.pdf']), '-dpdf', '-vector');
            else
                exportgraphics(fig_summary, fullfile(path_feature, ['summary_', feature_name, '.pdf']), 'ContentType', 'vector');
            end
        end
        % List all summary PDF files in the output path and its subdirectories.
        summary_files = dir(fullfile(path_out, '**', 'summary_*.pdf'));
        % Sort the summary PDF files by their folder names.
        [~, idx] = sort({summary_files.folder});
        summary_files = summary_files(idx);
        % Merge all the summary PDF files to a single PDF file.
        delete(fullfile(path_out, 'summary.pdf'));
        for i_file = 1:numel(summary_files)
            copyfile(fullfile(summary_files(i_file).folder, summary_files(i_file).name), path_out);
            mergePdfs(path_out, 'summary.pdf', path_out);
            delete(fullfile(path_out, summary_files(i_file).name));
        end
    end

    warning('on');

    % Close the figures.
    if n_rows > 0
        close(fig_summary);
        if ~profile_options.(ProfileOptionKey.SILENT.value) && profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
            fprintf('\nINFO: Summary stored in %s', path_out);
        end
    end

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('\nINFO: Finished the computation of the profiles under "%s" feature.\n', feature.name);
    end

    diary off;
end

% Following code is modified from the code provided by Benjamin Großmann (2024). Merge PDF-Documents
% (https://www.mathworks.com/matlabcentral/fileexchange/89127-merge-pdf-documents), MATLAB Central
% File Exchange. Retrieved November 12, 2024.
function mergePdfs(file_path, output_file_name, output_path)

    fileNames = dir(fullfile(file_path, '*.pdf'));
    fileNames = {fileNames.name};
    fileNames = cellfun(@(f) fullfile(file_path, f), fileNames, 'UniformOutput', false);

    memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
    merger = org.apache.pdfbox.multipdf.PDFMergerUtility;
    
    cellfun(@(f) merger.addSource(f), fileNames)
    
    merger.setDestinationFileName(fullfile(output_path, output_file_name));
    merger.mergeDocuments(memSet)
end

function integral = integrate(curve, profile_type, semilogx)
    % Compute the integral of the curve from different types of profiles.
    
    integral = 0;
    switch profile_type
        case 'perf'
            if semilogx
                curve(1, :) = log2(curve(1, :));
            end
            % The curve is a right-continuous step function.
            integral = integral + sum(diff(curve(1, :)) .* curve(2, 1:end-1));
        case 'data'
            if semilogx
                curve(1, :) = log2(1 + curve(1, :));
            end
            % The curve is a right-continuous step function.
            integral = integral + sum(diff(curve(1, :)) .* curve(2, 1:end-1));
        case 'log_ratio'
            integral = integral + sum(abs(curve(2, :)));
    end
end

function profile_scores = computeScores(curves, semilogx)
    % Compute the scores of the solvers for all the profiles.

    n_tols = size(curves, 2);
    n_solvers = size(curves{1}.hist.perf, 1);
    
    if n_solvers == 2
        profile_scores = ones(n_solvers, n_tols, 2, 3);
    else
        profile_scores = ones(n_solvers, n_tols, 2, 2);
    end

    for i_tol = 1:n_tols
        for i_solver = 1:n_solvers
            curve_hist_perf = curves{i_tol}.hist.perf{i_solver, end};
            curve_hist_data = curves{i_tol}.hist.data{i_solver, end};
            curve_out_perf = curves{i_tol}.out.perf{i_solver, end};
            curve_out_data = curves{i_tol}.out.data{i_solver, end};
            profile_scores(i_solver, i_tol, 1, 1) = integrate(curve_hist_perf, 'perf', semilogx);
            profile_scores(i_solver, i_tol, 1, 2) = integrate(curve_hist_data, 'data', semilogx);
            profile_scores(i_solver, i_tol, 2, 1) = integrate(curve_out_perf, 'perf', semilogx);
            profile_scores(i_solver, i_tol, 2, 2) = integrate(curve_out_data, 'data', semilogx);
            if n_solvers == 2
                curve_hist_log_ratio = curves{i_tol}.hist.log_ratio{i_solver};
                curve_out_log_ratio = curves{i_tol}.out.log_ratio{i_solver};
                if isempty(curve_hist_log_ratio)
                    curve_hist_log_ratio = [0; 0];
                end
                if isempty(curve_out_log_ratio)
                    curve_out_log_ratio = [0; 0];
                end
                profile_scores(i_solver, i_tol, 1, 3) = integrate(curve_hist_log_ratio, 'log_ratio', semilogx);
                profile_scores(i_solver, i_tol, 2, 3) = integrate(curve_out_log_ratio, 'log_ratio', semilogx);
            end
        end
    end

    % Normalize the profile_scores by dividing the maximum profile_scores with the same tolerance and types.
    max_scores = max(profile_scores, [], 1);
    max_scores(max_scores == 0) = 1;
    profile_scores = profile_scores ./ max_scores;
end