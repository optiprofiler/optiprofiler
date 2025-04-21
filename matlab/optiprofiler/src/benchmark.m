function [solver_scores, profile_scores, curves] = benchmark(varargin)
%BENCHMARK creates multiple profiles for benchmarking optimization solvers on a
%   set of problems with different features.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Signatures:
%
%   SOLVER_SCORES = BENCHMARK(SOLVERS) creates performance profiles, data
%   profiles, and log-ratio profiles for the given SOLVERS on the default
%   unconstrained problem set, returning SOLVER_SCORES based on the profiles.
%   SOLVERS is a cell array of function handles. We require SOLVERS to accept
%   specified inputs and return specified outputs. Details can be found in the
%   following 'Cautions' part.
%
%   SOLVER_SCORES = BENCHMARK(SOLVERS, FEATURE_NAME) creates profiles and data
%   profiles for the given SOLVERS on the default unconstrained problem set
%   with the specified feature FEATURE_NAME.
%
%   SOLVER_SCORES = BENCHMARK(SOLVERS, OPTIONS) creates profiles for the given
%   SOLVERS with options specified in the struct OPTIONS. See 'Options' part
%   for more details.
%
%   SOLVER_SCORES = BENCHMARK(OPTIONS) creates profiles with options specified
%   in the struct OPTIONS. Note that the struct OPTIONS should at least contain
%   the field `load` with the value 'latest' or a time stamp of an experiment
%   in the format of 'yyyyMMdd_HHmmss'. In this case, we will load the data
%   from the specified experiment and draw the profiles.
%
%   [SOLVER_SCORES, PROFILE_SCORES] = BENCHMARK(...) returns a 4D tensor
%   PROFILE_SCORES containing scores for all profiles. See `score_fun` in
%   'Options' part for more details.
%
%   [SOLVER_SCORES, PROFILE_SCORES, CURVES] = BENCHMARK(...) returns a cell
%   array CURVES containing the curves of all the profiles.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Options:
%
%   Options should be specified in a struct. The following are the available
%   fields of the struct:
%
%       1. options for profiles and plots:
%
%       - n_jobs: the number of parallel jobs to run the test. Default is 1.
%       - benchmark_id: the identifier of the test. It is used to create the
%         specific directory to store the results. Default is 'out'.
%       - feature_stamp: the stamp of the feature with the given options. It is
%         used to create the specific directory to store the results. Default
%         depends on features.
%       - range_type: the type of the uncertainty interval. For stochastic
%         features, we run several times of the experiments and get average
%         curves and uncertainty intervals. Default is 'minmax', meaning that
%         we takes the pointwise minimum and maximum of the curves.
%       - savepath: the path to store the results. Default is 'pwd', the
%         current working directory.
%       - max_tol_order: the maximum order of the tolerance. In any profile
%         (performance profiles, data profiles, and log-ratio profiles), we
%         need to set a group of 'tolerances' to define the convergence test of
%         the solvers. (Details can be found in the references.) We will set
%         the tolerances as `10^(-1:-1:-max_tol_order)`. Default is 10.
%       - max_eval_factor: the factor multiplied to each problem's dimension to
%         get the maximum number of evaluations for each problem. Default is
%         500.
%       - merit_fun: the merit function to measure the quality of a point using
%         the objective function value and the maximum constraint violation.
%         It should be a function handle as follows:
%               ``(fun_values, maxcv_values, maxcv_init) -> merit_values``,
%         where `fun_values` is history of the objective function values,
%         `maxcv_values` is history of the maximum constraint violation, and
%         `maxcv_init` is the initial maximum constraint violation. The size of
%         `fun_values` and `maxcv_values` is the same, and the size of
%         `maxcv_init` is the same as the second to last dimensions of
%         `fun_values`. The default merit function varphi(x) is defined by the
%         objective function f(x) and the maximum constraint violation v(x) as
%           varphi(x) = f(x)                        if v(x) <= v1
%           varphi(x) = f(x) + 1e5 * (v(x) - v1)    if v1 < v(x) <= v2
%           varphi(x) = Inf                         if v(x) > v2
%         where v1 = max(1e-5, v0) and v2 = min(0.01, 1e-10 * max(1, v0)),
%         and v0 is the initial maximum constraint violation.
%       - project_x0: whether to project the initial point to the feasible set.
%         Default is false.
%       - run_plain: whether to run an extra experiment with the 'plain'
%         feature. Default is false.
%       - score_only: whether to only calculate the scores of the solvers
%         without drawing the profiles and saving the data. Default is false.
%       - summarize_performance_profiles: whether to add all the performance
%         profiles to the summary PDF. Default is true.
%       - summarize_data_profiles: whether to add all the data profiles to the
%         summary PDF. Default is true.
%       - summarize_log_ratio_profiles: whether to add all the log-ratio
%         profiles to the summary PDF. Default is false.
%       - summarize_output_based_profiles: whether to add all the output-based
%         profiles of the selected profiles to the summary PDF. Default is
%         true.
%       - silent: whether to show the information of the progress. Default is
%         false.
%       - solver_verbose: the level of the verbosity of the solvers. `0` means
%         no verbosity, `1` means some verbosity, and `2` means full verbosity.
%         Default is 1.
%       - semilogx: whether to use the semilogx scale during plotting profiles
%         (performance profiles and data profiles). Default is true.
%       - normalized_scores: whether to normalize the scores of the solvers by
%         the maximum score of the solvers. Default is false.
%       - score_weight_fun: the weight function to calculate the scores of the
%         solvers in the performance and data profiles. It should be a function
%         handle representing a nonnegative function in R^+. Default is 1.
%       - score_fun: the scoring function to calculate the scores of the
%         solvers. It should be a function handle as follows:
%               ``profile_scores -> solver_scores``,
%         where `profile_scores` is a 4D tensor containing scores for all
%         profiles. The first dimension of `profile_scores` corresponds to the
%         index of the solver, the second corresponds to the index of tolerance
%         starting from 1, the third represents history-based or output-based
%         profiles, and the fourth represents performance profiles, data
%         profiles, or log-ratio profiles. The default scoring function takes
%         the average of the history-based performance profiles under all the
%         tolerances.
%       - load: loading the stored data from a completed experiment and draw
%         profiles. It can be either 'latest' or a time stamp of an experiment
%         in the format of 'yyyyMMdd_HHmmss'. No default.
%       - solvers_to_load: the indices of the solvers to load when the 'load'
%         option is provided. It can be a vector of different integers selected
%         from 1 to the total number of solvers of the loading experiment. At
%         least two indices should be provided. Default is all the solvers.
%       - line_colors: the colors of the lines in the plots. It can be a cell
%         array of short names of colors ('r', 'g', 'b', 'c', 'm', 'y', 'k') or
%         a matrix with each row being a RGB triplet. Default line colors are
%         those in the palettename named "gem" (see MATLAB documentation for
%         'colororder'). Note that if the number of solvers is greater than the
%         number of colors, we will cycle through the colors.
%       - line_styles: the styles of the lines in the plots. It can be a cell
%         array of chars that are the combinations of line styles ('-', '-.',
%         ':', '--') and markers ('none', 'o', '+', '*', '.', 'x', 's', 'd',
%         '^', 'v', '>', '<', 'p', 'h'). Default line style order is {'-',
%         '-.', ':', '--'}. Note that if the number of solvers is greater than
%         the number of line styles, we will cycle through the styles.
%       - line_widths: the widths of the lines in the plots. It should be a
%         positive scalar or a vector. Default is 1.5. Note that if the number
%         of solvers is greater than the number of line widths, we will cycle
%         through the widths.
%       - bar_colors: two different colors for the bars of two solvers in the
%         log-ratio profiles. It can be a cell array of short names of colors
%         ('r', 'g', 'b', 'c', 'm', 'y', 'k') or a 2-by-3 matrix with each row
%         being a RGB triplet. Default is set to the first two colors in the
%         'line_colors' option.
%
%       2. options for features:
%
%       - feature_name: the name of the feature. The available features are
%         'plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted',
%         'linearly_transformed', 'random_nan', 'unrelaxable_constraints',
%         'nonquantifiable_constraints', 'quantized', and 'custom'. Default is
%         'plain'.
%       - n_runs: the number of runs of the experiments under the given
%         feature. Default is 10 for stochastic features and 1 for
%         deterministic features.
%       - distribution: the distribution of perturbation in 'perturbed_x0'
%         feature or noise in 'noisy' feature. It should be either a string
%         (or char), or a function handle
%               ``(random_stream, dimension) -> random vector``,
%         accepting a random_stream and the dimension of a problem and
%         returning a random vector with the given dimension. In 'perturbed_x0'
%         case, the char should be either 'spherical' or 'gaussian' (default is
%         'spherical'). In 'noisy' case, the char should be either 'gaussian'
%         or 'uniform' (default is 'gaussian').
%       - perturbation_level: the magnitude of the perturbation to the initial
%         guess in the 'perturbed_x0' feature. Default is 1e-3.
%       - noise_level: the magnitude of the noise in the 'noisy' feature.
%         Default is 1e-3.
%       - noise_type: the type of the noise in the 'noisy' features. It should
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
%         be 2 ^ (condition_factor * n / 2), where `n` is the dimension of the
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
%         is 1e-3.
%       - mesh_type: the type of the mesh in the 'quantized' feature. It should
%         be either 'absolute' or 'relative'. Default is 'absolute'.
%       - ground_truth: whether the featured problem is the ground truth or not
%         in the 'quantized' feature. Default is true.
%       - mod_x0: the modifier function to modify the inital guess in the 
%         'custom' feature. It should be a function handle as follows:
%               ``(random_stream, problem) -> modified_x0``,
%         where `problem` is an instance of the class Problem, and
%         `modified_x0` is the modified initial guess. No default.
%       - mod_affine: the modifier function to generate the affine
%         transformation applied to the variables in the 'custom' feature. It
%         should be a function handle as follows:
%               ``(random_stream, problem) -> (A, b, inv)``,
%         where `problem` is an instance of the class Problem, `A` is the
%         matrix of the affine transformation, `b` is the vector of the affine
%         transformation, and `inv` is the inverse of matrix `A`. No default.
%       - mod_bounds: the modifier function to modify the bound constraints in
%         the 'custom' feature. It should be a function handle as follows:
%               ``(random_stream, problem) -> (modified_xl, modified_xu)``,
%         where `problem` is an instance of the class Problem, `modified_xl` is
%         the modified lower bound, and `modified_xu` is the modified upper
%         bound. No default.
%       - mod_linear_ub: the modifier function to modify the linear inequality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(random_stream, problem) -> (modified_aub, modified_bub)``,
%         where `problem` is an instance of the class Problem, `modified_aub`
%         is the modified matrix of the linear inequality constraints, and
%         `modified_bub` is the modified vector of the linear inequality
%         constraints. No default.
%       - mod_linear_eq: the modifier function to modify the linear equality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(random_stream, problem) -> (modified_aeq, modified_beq)``,
%         where `problem` is an instance of the class Problem, `modified_aeq`
%         is the modified matrix of the linear equality constraints, and
%         `modified_beq` is the modified vector of the linear equality
%         constraints. No default.
%       - mod_fun: the modifier function to modify the objective function in
%         the 'custom' feature. It should be a function handle as follows:
%               ``(x, random_stream, problem) -> modified_fun``,
%         where `x` is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_fun` is the modified objective function
%         value. No default.
%       - mod_cub: the modifier function to modify the nonlinear inequality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(x, random_stream, problem) -> modified_cub``,
%         where x is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_cub` is the modified vector of the
%         nonlinear inequality constraints. No default.
%       - mod_ceq: the modifier function to modify the nonlinear equality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(x, random_stream, problem) -> modified_ceq``,
%         where x is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_ceq` is the modified vector of the
%         nonlinear equality constraints. No default.
%
%       3. options for CUTEst:
%
%       Note that the CUTEst we used is the MATLAB codes from a GitHub
%       repository called 'S2MPJ', created by Professor Serge Gratton and
%       Professor Philippe L. Toint. More details can be found in the following
%       website.
%           https://github.com/GrattonToint/S2MPJ
%
%       - ptype: the type of the problems to be selected. It should be a string
%         or char consisting of any combination of 'u' (unconstrained), 'b'
%         (bound constrained), 'l' (linearly constrained), and 'n' (nonlinearly
%         constrained), such as 'b', 'ul', 'ubn'. Default is 'u'.
%       - mindim: the minimum dimension of the problems to be selected. Default
%         is 1.
%       - maxdim: the maximum dimension of the problems to be selected. Default
%         is mindim + 1.
%       - mincon: the minimum number of linear and nonlinear constraints of the
%         problems to be selected. Default is 0.
%       - maxcon: the maximum number of linear and nonlinear constraints of the
%         problems to be selected. Default is mincon + 10.
%       - excludelist: the list of problems to be excluded. Default is not to
%         exclude any problem.
%
%       Note that if the `load` option is provided, we will use above options
%       to select the problems and then take an intersection with the problems
%       in the loading experiment.
%
%       4. other options:
%
%       - solver_names: the names of the solvers. Default is the names of the
%         function handles in `solvers`.
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
%               ``(problem_name) -> problem``,
%         where `problem_name` is the name of the problem, and `problem` is an
%         instance of the class Problem. Default is not to load any custom
%         problem.
%       - custom_problem_names: the names of the custom problems to be
%         selected. Default is not to select any custom problem.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cautions:
%
%   1. Each solver in SOLVERS should accept the following signature(s):
%
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
%
%   2. The log-ratio profiles are available only when there are exactly two 
%      solvers.
%
%   For more information of performance and data profiles, see [1]_, [2]_,
%   [5]_. For that of log-ratio profiles, see [4]_, [6]_. For that of S2MPJ,
%   see [3]_.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   References:
%
%   .. [1] E. D. Dolan and J. J. Moré. Benchmarking optimization software with
%          performance profiles. *Math. Program.*, 91(2):201–213, 2002.
%          doi:10.1007/s101070100263
%          <https://doi.org/10.1007/s101070100263>.
%   .. [2] N. Gould and J. Scott. A note on performance profiles for
%          benchmarking software. *ACM Trans. Math. Software*, 43(2):15:1–5,
%          2016. doi:10.1145/2950048 <https://doi.org/10.1145/2950048>.
%   .. [3] S. Gratton and Ph. L. Toint. S2MPJ and CUTEst optimization problems
%          for Matlab, Python and Julia. arXiv:2407.07812, 2024.
%   .. [4] J. L. Morales. A numerical study of limited memory BFGS methods.
%          *Appl. Math. Lett.*, 15(4):481–487, 2002.
%          doi:10.1016/S0893-9659(01)00162-8
%          <https://doi.org/10.1016/S0893-9659(01)00162-8>.
%   .. [5] J. J. Moré and S. M. Wild. Benchmarking derivative-free optimization
%          algorithms. *SIAM J. Optim.*, 20(1):172–191, 2009.
%          doi:10.1137/080724083 <https://doi.org/10.1137/080724083>.
%   .. [6] H.-J. M. Shi, M. Q. Xuan, F. Oztoprak, and J. Nocedal. On the
%          numerical performance of finite-difference-based methods for
%          derivative-free optimization. *Optim. Methods Softw.*,
%          38(2):289–311, 2023. doi:10.1080/10556788.2022.2121832
%          <https://doi.org/10.1080/10556788.2022.2121832>.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Authors:
%               Cunxin HUANG (cun-xin.huang@connect.polyu.hk)
%               Tom M. RAGONNEAU (t.ragonneau@gmail.com)
%               Zaikun ZHANG (zhangzaikun@mail.sysu.edu.cn)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the input arguments. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 0
        error("MATLAB:benchmark:solverMustBeProvided", "A cell of function handles (callable solvers) or a struct of options must be provided.");
    elseif nargin == 1
        if isstruct(varargin{1})
            % When input contains one argument and the first argument is a struct, we assume the
            % user chooses benchmark(options).
            solvers = {};
            feature_name = 'plain';
            options = varargin{1};
            options_store = options;
            if ~isfield(options, ProfileOptionKey.LOAD.value) || isempty(options.(ProfileOptionKey.LOAD.value))
                error("MATLAB:benchmark:LoadFieldNotProvided", "The field 'load' of options should be provided when the first argument is a struct.");
            end
        else
            solvers = varargin{1};
            feature_name = 'plain';
            options = struct();
        end
    elseif nargin == 2
        if ischarstr(varargin{2})
            % When input contains two arguments and the second argument is a char or cell of char,
            % we assume the user chooses benchmark(solvers, feature_name).
            solvers = varargin{1};
            feature_name = varargin{2};
            options = struct();
        elseif isstruct(varargin{2})
            % When input contains two arguments and the second argument is a struct, we assume the
            % user chooses benchmark(solvers, options).
            solvers = varargin{1};
            options = varargin{2};
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Process the 'load' option if it is provided. %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isfield(options, ProfileOptionKey.LOAD.value) && ~isempty(options.(ProfileOptionKey.LOAD.value))
        % Check the validity of the 'load' field.
        if ~ischarstr(options.(ProfileOptionKey.LOAD.value))
            error("MATLAB:checkValidityProfileOptions:loadNotValid", "The field 'load' of options should be a char or a string.");
        end
        options.(ProfileOptionKey.LOAD.value) = char(options.(ProfileOptionKey.LOAD.value));
        % Check whether it is 'latest' or a string (or char) in the format of 'yyyyMMdd_HHmmss'.
        time_stamp_pattern = '^\d{4}(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})$';
        if ~strcmp(options.(ProfileOptionKey.LOAD.value), 'latest') && isempty(regexp(options.(ProfileOptionKey.LOAD.value), time_stamp_pattern, 'once'))
            error("MATLAB:checkValidityProfileOptions:loadNotValid", "The field 'load' of options should be either 'latest' or a time stamp in the format of 'yyyyMMdd_HHmmss'.");
        end

        % Set the path to search for the data to load.
        if isfield(options, ProfileOptionKey.BENCHMARK_ID.value)
            search_path = fullfile(pwd, options.(ProfileOptionKey.BENCHMARK_ID.value));
        else
            search_path = pwd;
            options.(ProfileOptionKey.BENCHMARK_ID.value) = '.';
        end

        % Find the path of the data to load.
        time_stamp_files = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
        if strcmp(options.(ProfileOptionKey.LOAD.value), 'latest')
            % Try to find all the txt files named by time_stamp in the search_path directory and find the latest one.
            % Note that we limit the search to 5 levels of subdirectories.
            time_stamp_files = search_in_dir(search_path, 'time_stamp_*.txt', 5, 0, time_stamp_files);
            if isempty(time_stamp_files)
                error("MATLAB:benchmark:NoTimeStamps", "Failed to load data since no time_stamp files are found in the directory '%s'. Note that the search is limited to 5 levels of subdirectories.", search_path);
            end
            time_stamps = arrayfun(@(f) datetime(f.name(12:end-4), 'InputFormat', 'yyyyMMdd_HHmmss'), time_stamp_files);
            [~, indexes] = sort(time_stamps, 'descend');
            % Get the latest time_stamp file. If there are multiple files with the same time_stamp, we take the first one.
            latest_idx = indexes(1);
            latest_time_stamp_file = time_stamp_files(latest_idx);
            path_data = latest_time_stamp_file.folder;
        else
            % Same as above, but we only search for the specific time_stamp file.
            pattern = ['time_stamp_', options.(ProfileOptionKey.LOAD.value), '.txt'];
            time_stamp_file = search_in_dir(search_path, pattern, 5, 0, time_stamp_files);
            if isempty(time_stamp_file)
                error("MATLAB:benchmark:NoTimeStamps", "Failed to load data since no time_stamp named '%s' is found in the directory '%s'. Note that the search is limited to 5 levels of subdirectories.", ['time_stamp_', options.(ProfileOptionKey.LOAD.value), '.txt'], search_path);
            end
            path_data = time_stamp_file.folder;
        end

        % Load data from the 'data.mat' file in the path_data directory.
        if ~isfield(options, ProfileOptionKey.FEATURE_STAMP.value)
            warning('off');
            load(fullfile(path_data, 'data.mat'), 'feature_stamp');
            warning('on');
            if ~exist('feature_stamp', 'var')
                error("MATLAB:benchmark:NoFeatureStampDataMatFile", "Failed to load the variable 'feature_stamp' from the 'data.mat' file in the directory '%s' and the 'feature_stamp' field is not provided in the options.", path_data);
            end
        end
        warning('off');
        load(fullfile(path_data, 'data.mat'), 'n_evals', 'problem_names', 'problem_types', 'problem_dims', 'problem_cons', 'computation_times', 'solvers_successes', 'merit_histories', 'merit_out', 'merit_init');
        warning('on');
        if ~exist('n_evals', 'var') || ~exist('problem_names', 'var') || ~exist('problem_types', 'var') || ~exist('problem_dims', 'var') || ~exist('problem_cons', 'var') || ~exist('computation_times', 'var') || ~exist('solvers_successes', 'var') || ~exist('merit_histories', 'var') || ~exist('merit_out', 'var') || ~exist('merit_init', 'var')
            error("MATLAB:benchmark:NoDataMatFile", "Failed to load computation results from the 'data.mat' file in the directory '%s'.", path_data);
        end

        % Check whether the user provides the 'solvers_to_load' field.
        n_solvers_loaded = size(merit_out, 2);
        if isfield(options, ProfileOptionKey.SOLVERS_TO_LOAD.value)
            % Check the validity of the 'solvers_to_load' field. It should be a vector of different
            % integers selected from 1 to the total number of solvers in the loaded data.
            try
                solvers_to_load = unique(options.(ProfileOptionKey.SOLVERS_TO_LOAD.value));
            catch
            end
            if ~isintegervector(solvers_to_load) || any(solvers_to_load < 1) || any(solvers_to_load > n_solvers_loaded) || numel(solvers_to_load) < 2
                error("MATLAB:benchmark:solvers_to_loadNotValid", "The field 'solvers_to_load' of options should be a vector of different integers selected from 1 to the total number of solvers in the loaded data, and at least two indices should be provided.");
            end
        else
            solvers_to_load = 1:n_solvers_loaded;
        end

        % Truncate the loaded data according to the 'solvers_to_load' field.
        if ~isfield(options, OtherOptionKey.SOLVER_NAMES.value)
            warning('off');
            load(fullfile(path_data, 'data.mat'), 'solver_names');
            warning('on');
            if ~exist('solver_names', 'var')
                error("MATLAB:benchmark:NoSolverNamesDataMatFile", "Failed to load the variable 'solver_names' from the 'data.mat' file in the directory '%s' and the 'solver_names' field is not provided in the options.", path_data);
            end
            solver_names = solver_names(solvers_to_load);
        elseif numel(options.(OtherOptionKey.SOLVER_NAMES.value)) ~= numel(solvers_to_load)
            error("MATLAB:benchmark:solver_namesNotMatch", "The number of elements in the 'solver_names' field of options should be the same as the number of solvers to load in the 'solvers_to_load' field.");
        end
        n_solvers = numel(solvers_to_load);
        n_evals = n_evals(:, solvers_to_load, :);
        computation_times = computation_times(:, solvers_to_load, :);
        solvers_successes = solvers_successes(:, solvers_to_load, :);
        merit_histories = merit_histories(:, solvers_to_load, :, :);
        merit_out = merit_out(:, solvers_to_load, :);

        % Check whether the user provides options for CUTEst. If so, we will use them to filter the
        % problems in the loaded data.
        options = checkValidityCutestOptions(options);
        p_to_load = true(size(problem_names));
        if isfield(options, CutestOptionKey.PTYPE.value)
            % Judge whether the type of the problems in the loaded data satisfies the options.
            for i = 1:numel(problem_types)
                if ~ismember(problem_types{i}, options.(CutestOptionKey.PTYPE.value))
                    p_to_load(i) = false;
                end
            end
        end
        if isfield(options, CutestOptionKey.MINDIM.value)
            % Judge whether the dimension of the problems in the loaded data satisfies the options.
            for i = 1:numel(problem_dims)
                if problem_dims(i) < options.(CutestOptionKey.MINDIM.value)
                    p_to_load(i) = false;
                end
            end
        end
        if isfield(options, CutestOptionKey.MAXDIM.value)
            % Judge whether the dimension of the problems in the loaded data satisfies the options.
            for i = 1:numel(problem_dims)
                if problem_dims(i) > options.(CutestOptionKey.MAXDIM.value)
                    p_to_load(i) = false;
                end
            end
        end
        if isfield(options, CutestOptionKey.MINCON.value)
            % Judge whether the number of constraints of the problems in the loaded data satisfies
            % the options.
            for i = 1:numel(problem_cons)
                if problem_cons(i) < options.(CutestOptionKey.MINCON.value)
                    p_to_load(i) = false;
                end
            end
        end
        if isfield(options, CutestOptionKey.MAXCON.value)
            % Judge whether the number of constraints of the problems in the loaded data satisfies
            % the options.
            for i = 1:numel(problem_cons)
                if problem_cons(i) > options.(CutestOptionKey.MAXCON.value)
                    p_to_load(i) = false;
                end
            end
        end
        if isfield(options, CutestOptionKey.EXCLUDELIST.value)
            % Judge whether the problems in the loaded data are in the exclude list.
            for i = 1:numel(problem_names)
                if ismember(problem_names{i}, options.(CutestOptionKey.EXCLUDELIST.value))
                    p_to_load(i) = false;
                end
            end
        end
        
        if any(p_to_load)
            % Truncate the loaded data according to the selected problems.
            problem_names = problem_names(p_to_load);
            n_evals = n_evals(p_to_load, :, :);
            computation_times = computation_times(p_to_load, :, :);
            solvers_successes = solvers_successes(p_to_load, :, :);
            merit_histories = merit_histories(p_to_load, :, :, :);
            merit_out = merit_out(p_to_load, :, :);
        else
            error("MATLAB:benchmark:NoProblemsToLoad", "Failed to load data since no problems are selected by the options.");
        end

        merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Process `solvers` and `feature_name` %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process the solvers.
    if ~isfield(options, ProfileOptionKey.LOAD.value) || isempty(options.(ProfileOptionKey.LOAD.value))
        if ~iscell(solvers) || ~all(cellfun(@(s) isa(s, 'function_handle'), solvers))
            error("MATLAB:benchmark:solversWrongType", "The first argument for `benchmark` must be a cell array of function handles.");
        end
        if numel(solvers) < 2
            error("MATLAB:benchmark:solversAtLeastTwo", "The first argument for `benchmark` must be a cell array of at least two function handles since we need to compare at least two solvers.");
        end
    end

    % Process the feature_name.
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

    if ~exist('n_solvers', 'var')
        n_solvers = numel(solvers);
    end
    solver_scores = zeros(n_solvers, 1);
    profile_scores = [];
    curves = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Use cutest_options to select problems. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = loadCutestNames(cutest_options, other_options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Set the default options for plotting. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Define the directory and load the results. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define the directory to store or load the results.
    path_out = fullfile(profile_options.(ProfileOptionKey.SAVEPATH.value), profile_options.(ProfileOptionKey.BENCHMARK_ID.value));

    % Record the time stamp of the current experiment.
    time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
    time_stamp = char(time);

    % Set the default feature stamp if it does not exist.
    if ~exist('feature_stamp', 'var')
        feature_stamp = profile_options.(ProfileOptionKey.FEATURE_STAMP.value);
    end

    path_feature = fullfile(path_out, [feature_stamp, '_', time_stamp]);
    path_log = fullfile(path_feature, 'test_log');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Create the directory to store the results. %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        if ~exist(path_feature, 'dir')
            mkdir(path_feature);
        else
            rmdir(path_feature, 's');
            mkdir(path_feature);
        end
    end

    if profile_options.(ProfileOptionKey.SCORE_ONLY.value) || ~isempty(profile_options.(ProfileOptionKey.LOAD.value))
        % If 'load' is not empty, we will directly load the results and do not need to compute
        % again. In this case, we do not need to create the directory to store the hist plots.
        path_hist_plots = '';
    else
        path_hist_plots = fullfile(path_feature, 'history_plots');
        if ~exist(path_hist_plots, 'dir')
            mkdir(path_hist_plots);
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Create the directory to store options and log files.
        if ~exist(path_log, 'dir')
            mkdir(path_log);
        else
            rmdir(path_log, 's');
            mkdir(path_log);
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Create a README.txt file to explain the content of the folder `path_log`.
        path_readme_log = fullfile(path_log, 'README.txt');
        try
            fid = fopen(path_readme_log, 'w');
            fprintf(fid, "This folder contains following files and directories.\n\n");
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to create the README.txt file for folder '%s'.\n", path_log);
            end
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Create a txt file named by time_stamp to record the time_stamp.
        path_time_stamp = fullfile(path_log, ['time_stamp_', time_stamp, '.txt']);
        try
            fid = fopen(path_time_stamp, 'w');
            fprintf(fid, time_stamp);
            fclose(fid);
            try
                fid = fopen(path_readme_log, 'a');
                fprintf(fid, "'time_stamp_%s.txt': file, recording the time stamp of the current experiment.\n", time_stamp);
                fclose(fid);
            catch
            end
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to create the time_stamp file for folder '%s'.\n", path_log);
            end
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value) && ~isfield(other_options, OtherOptionKey.PROBLEM.value)
        path_figs = fullfile(path_log, 'figs');
        mkdir(path_figs);
        try
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'figs': folder, containing all the FIG files of the profiles.\n");
            fclose(fid);
        catch
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        try
            if exist('options_store', 'var')
                save(fullfile(path_log, 'options_store.mat'), 'options_store');
                try
                    fid = fopen(path_readme_log, 'a');
                    fprintf(fid, "'options_store.mat': file, storing the options of the current experiment.\n");
                    fclose(fid);
                catch
                end
            end
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to save the `options` of the current experiment.\n");
            end
        end
        log_file = fullfile(path_log, 'log.txt');
        diary(log_file);
        try
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'log.txt': file, the log file of the current experiment, recording print information from MATLAB command window.\n");
            fclose(fid);
        catch
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Create a README.txt file to explain the content of the folder `path_feature`.
        path_readme_feature = fullfile(path_feature, 'README.txt');
        try
            fid = fopen(path_readme_feature, 'w');
            fprintf(fid, "This folder contains following files and directories.\n\n");
            if ~isempty(path_hist_plots)
                fprintf(fid, "'history_plots': folder, containing all the history plots for each problem.\n");
                fprintf(fid, "'history_plots_summary.pdf': file, the summary PDF of history plots for all problems.\n");
            end
            fprintf(fid, "'test_log': folder, containing log files and other useful experimental data.\n");
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to create the README.txt file for folder '%s'.\n", path_feature);
            end
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % We try to copy the script or function that calls the benchmark function to the log directory.
        try
            calling_script = dbstack(1, '-completenames');
            if ~isempty(calling_script)
                copyfile(calling_script.file, path_log);
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fprintf("INFO: The script or function that calls `benchmark` function is copied to: %s.\n\n", path_log);
                end
                try
                    fid = fopen(path_readme_log, 'a');
                    fprintf(fid, "'%s': file, the script or function that calls `benchmark` function.\n", calling_script.file);
                    fclose(fid);
                catch
                end
            end
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to copy the script or function that calls `benchmark` function to the log directory.\n\n");
            end
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value) && ~isfield(other_options, OtherOptionKey.PROBLEM.value)
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
        if n_solvers == 2
            if ~exist(path_log_ratio_hist, 'dir')
                mkdir(path_log_ratio_hist);
            end
            if ~exist(path_log_ratio_out, 'dir')
                mkdir(path_log_ratio_out);
            end
        end

        try
            fid = fopen(path_readme_feature, 'a');
            fprintf(fid, "'detailed_profiles': folder, containing all the high-quality single profiles.\n");
            fclose(fid);
        catch
        end

        pdf_perf_hist_summary = fullfile(path_feature, 'perf_hist.pdf');
        pdf_perf_out_summary = fullfile(path_feature, 'perf_out.pdf');
        pdf_data_hist_summary = fullfile(path_feature, 'data_hist.pdf');
        pdf_data_out_summary = fullfile(path_feature, 'data_out.pdf');
        pdf_log_ratio_hist_summary = fullfile(path_feature, 'log-ratio_hist.pdf');
        pdf_log_ratio_out_summary = fullfile(path_feature, 'log-ratio_out.pdf');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Solve all the problems when 'load' option is not provided. %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % If `n_runs` is not specified, the feature is deterministic, and at least one solver is
    % randomized, then we set `n_runs` to 5.
    if ~isfield(feature_options, FeatureOptionKey.N_RUNS.value) && ~feature.is_stochastic && any(other_options.(OtherOptionKey.SOLVER_ISRAND.value))
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf("INFO: We set `n_runs` to 5 since the feature is deterministic and at least one solver is randomized and `n_runs` is not specified.\n\n");
        end
        feature.options.(FeatureOptionKey.N_RUNS.value) = 5;
    end

    if isempty(profile_options.(ProfileOptionKey.LOAD.value))
        % If 'load' is not specified, we solve all the problems.
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('INFO: Starting the computation of the "%s" profiles.\n', feature.name);
        end
        [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_evals, problem_names, problem_types, problem_dims, problem_cons,computation_times, solvers_successes] = solveAllProblems(solvers, feature, profile_options, other_options, true, path_hist_plots);
        merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
        try
            merit_histories = merit_fun(fun_histories, maxcv_histories, maxcv_init);
            merit_out = merit_fun(fun_out, maxcv_out, maxcv_init);
            merit_init = merit_fun(fun_init, maxcv_init, maxcv_init);
        catch
            error("MATLAB:benchmark:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
        end
        merit_min = min(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');

        % If there are no problems solved, skip the rest of the code, print a message, and return.
        if isempty(problem_names) || all(~solvers_successes(:))
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf('INFO: No problems were solved for the "%s" feature.\n', feature.name);
            end
            if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
                diary off;
            end
            return;
        end

        % If a specific problem is provided to `other_options`, we only solve this problem and generate
        % the history plots for it.
        if isfield(other_options, OtherOptionKey.PROBLEM.value)
            if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
                % We move the history plots to the feature directory.
                try
                    movefile(fullfile(path_hist_plots, '*'), path_feature);
                    rmdir(path_hist_plots, 's');
                    if ~profile_options.(ProfileOptionKey.SILENT.value)
                        fprintf('\nINFO: Detailed results stored in: \n%s\n\n', path_feature);
                    end
                catch
                end
            end

            % Since we will not compute the profiles, we set `solver_scores` to be the relative
            % decrease in the objective function value.
            solver_merit_mins = squeeze(min(min(merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'));
            solver_scores = (merit_init - solver_merit_mins) ./ max(merit_init - merit_min, eps);
            solver_scores = solver_scores';

            if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
                diary off;
            end
            return;
        end

        % Determine the least merit value for each problem.
        if feature.run_plain && profile_options.(ProfileOptionKey.RUN_PLAIN.value)
            feature_plain = Feature(FeatureName.PLAIN.value);
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf('\nINFO: Starting the computation of the profiles under "plain" feature.\n');
            end
            [fun_histories_plain, maxcv_histories_plain, ~, ~, ~, ~, ~, problem_names_plain, ~, ~, ~, computation_times_plain] = solveAllProblems(solvers, feature_plain, profile_options, other_options, false, {});
            merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
            try
                merit_histories_plain = merit_fun(fun_histories_plain, maxcv_histories_plain, maxcv_init);
            catch
                error("MATLAB:benchmark:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
            end
            merit_min_plain = min(min(min(merit_histories_plain, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
            computation_times = cat(3, computation_times, NaN(size(computation_times_plain)));
            
            for i_problem = 1:numel(problem_names)
                idx = find(strcmp(problem_names{i_problem}, problem_names_plain), 1);
                % Redefine the merit_min for the problems that are solved under the plain feature.
                % Note that min(x, NaN) = x.
                merit_min(i_problem) = min(merit_min(i_problem), merit_min_plain(idx), 'omitnan');
                computation_times(i_problem, :, size(computation_times, 3)) = computation_times_plain(idx, :, :);
            end
        end

        if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
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
                fprintf('\nINFO: Detailed results stored in: \n%s\n\n', path_feature);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Store the problem names. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Pick out unsolved problems.
    unsolved_problems = [];
    for i_problem = 1:numel(problem_names)
        if all(~solvers_successes(i_problem, :, :))
            unsolved_problems = [unsolved_problems, problem_names(i_problem)];
        end
    end

    % Calculate the computation times for each problem.
    time_processes = zeros(numel(problem_names), 1);
    for i_problem = 1:numel(problem_names)
        time_process = computation_times(i_problem, :, :);
        time_processes(i_problem) = sum(time_process(:), 'omitnan');
    end

    % Store the names of the problems.
    path_txt = fullfile(path_log, 'problems.txt');
    [~, idx] = sort(lower(problem_names));
    sorted_problem_names = problem_names(idx);
    sorted_time_processes = time_processes(idx);
    max_name_length = max(cellfun(@length, sorted_problem_names));
    sorted_time_processes = num2cell(sorted_time_processes);
    max_time_length = max(cellfun(@(x) length(sprintf('%.2f', x)), sorted_time_processes));

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        try
            fid = fopen(path_txt, 'w');
            for i = 1:length(sorted_problem_names)
                if ismember(sorted_problem_names{i}, unsolved_problems)
                    continue;
                end
                count = fprintf(fid, "%-*s      %*s\n", max_name_length, sorted_problem_names{i}, max_time_length, sprintf('%.2f seconds', sorted_time_processes{i}));
                if count < 0
                    if ~profile_options.(ProfileOptionKey.SILENT.value)
                        fprintf("INFO: Failed to record data for %s.", sorted_problem_names{i});
                    end
                end
            end
            if ~isempty(unsolved_problems)
                fprintf(fid, "\n");
                fprintf(fid, "Unsolved problems (that all the solvers failed to return a solution in all runs):\n");
                for i = 1:length(unsolved_problems)
                    count = fprintf(fid, "%s\n", unsolved_problems{i});
                    if count < 0
                        if ~profile_options.(ProfileOptionKey.SILENT.value)
                            fprintf("INFO: Failed to record data for %s.", unsolved_problems{i});
                        end
                    end
                end
            end
            fclose(fid);
            try
                fid = fopen(path_readme_log, 'a');
                fprintf(fid, "'problems.txt': file, storing the names of the solved and unsolved problems, the time spent on solving each problem, and the names of the problems that all the solvers failed to meet the convergence test for every tolerance and run.\n");
                fclose(fid);
            catch
            end
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Error occurred when writing the problem names to %s.\n", path_txt);
            end
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
        T_feature_stamp = strrep(feature_stamp, '_', '\_');
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
    curves = cell(1, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('\n');
    end

    % Find the problems that all the solvers failed to meet the convergence test for every tolerance.
    solvers_all_diverge_hist = false(n_problems, n_runs, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));
    solvers_all_diverge_out = false(n_problems, n_runs, profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value));

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
                        threshold = max(tolerance * merit_init(i_problem) + (1 - tolerance) * merit_min(i_problem), merit_min(i_problem));
                    else
                        threshold = -Inf;
                    end
                    if min(merit_histories(i_problem, i_solver, i_run, :), [], 'omitnan') <= threshold
                        work_hist(i_problem, i_solver, i_run) = find(merit_histories(i_problem, i_solver, i_run, :) <= threshold, 1, 'first');
                    end
                    if merit_out(i_problem, i_solver, i_run) <= threshold
                        work_out(i_problem, i_solver, i_run) = n_evals(i_problem, i_solver, i_run);
                    end
                end
            end
        end

        for i_problem = 1:n_problems
            for i_run = 1:n_runs
                solvers_all_diverge_hist(i_problem, i_run, i_tol) = all(isnan(work_hist(i_problem, :, i_run)));
                solvers_all_diverge_out(i_problem, i_run, i_tol) = all(isnan(work_out(i_problem, :, i_run)));
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

        if ~exist('solver_names', 'var')
            solver_names = cellfun(@(s) strrep(s, '_', '\_'), other_options.(OtherOptionKey.SOLVER_NAMES.value), 'UniformOutput', false);
        end

        [fig_perf_hist, fig_data_hist, fig_log_ratio_hist, curves{i_tol}.hist] = drawProfiles(work_hist, problem_dims, solver_names, tolerance_label, cell_axs_summary_hist, true, is_perf, is_data, is_log_ratio, profile_options, curves{i_tol}.hist);
        [fig_perf_out, fig_data_out, fig_log_ratio_out, curves{i_tol}.out] = drawProfiles(work_out, problem_dims, solver_names, tolerance_label, cell_axs_summary_out, is_output_based, is_perf, is_data, is_log_ratio, profile_options, curves{i_tol}.out);

        if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
            pdf_perf_hist = fullfile(path_perf_hist, ['perf_hist_', int2str(i_tol), '.pdf']);
            figure_perf_hist = fullfile(path_figs, ['perf_hist_', int2str(i_tol), '.fig']);
            exportgraphics(fig_perf_hist, pdf_perf_hist, 'ContentType', 'vector');
            savefig(fig_perf_hist, figure_perf_hist);
            pdf_perf_out = fullfile(path_perf_out, ['perf_out_', int2str(i_tol), '.pdf']);
            figure_perf_out = fullfile(path_figs, ['perf_out_', int2str(i_tol), '.fig']);
            exportgraphics(fig_perf_out, pdf_perf_out, 'ContentType', 'vector');
            savefig(fig_perf_out, figure_perf_out);
            pdf_data_hist = fullfile(path_data_hist, ['data_hist_', int2str(i_tol), '.pdf']);
            figure_data_hist = fullfile(path_figs, ['data_hist_', int2str(i_tol), '.fig']);
            exportgraphics(fig_data_hist, pdf_data_hist, 'ContentType', 'vector');
            savefig(fig_data_hist, figure_data_hist);
            pdf_data_out = fullfile(path_data_out, ['data_out_', int2str(i_tol), '.pdf']);
            figure_data_out = fullfile(path_figs, ['data_out_', int2str(i_tol), '.fig']);
            exportgraphics(fig_data_out, pdf_data_out, 'ContentType', 'vector');
            savefig(fig_data_out, figure_data_out);
            if n_solvers <= 2
                pdf_log_ratio_hist = fullfile(path_log_ratio_hist, ['log-ratio_hist_', int2str(i_tol), '.pdf']);
                figure_log_ratio_hist = fullfile(path_figs, ['log-ratio_hist_', int2str(i_tol), '.fig']);
                exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist, 'ContentType', 'vector');
                savefig(fig_log_ratio_hist, figure_log_ratio_hist);
                pdf_log_ratio_out = fullfile(path_log_ratio_out, ['log-ratio_out_', int2str(i_tol), '.pdf']);
                figure_log_ratio_out = fullfile(path_figs, ['log-ratio_out_', int2str(i_tol), '.fig']);
                exportgraphics(fig_log_ratio_out, pdf_log_ratio_out, 'ContentType', 'vector');
                savefig(fig_log_ratio_out, figure_log_ratio_out);
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
                if n_solvers == 2
                    exportgraphics(fig_log_ratio_hist, pdf_log_ratio_hist_summary, 'ContentType', 'vector', 'Append', true);
                end
                if ~isempty(fig_log_ratio_out)
                    exportgraphics(fig_log_ratio_out, pdf_log_ratio_out_summary, 'ContentType', 'vector', 'Append', true);
                end
            end

            if i_tol == 1
                try
                    fid = fopen(path_readme_feature, 'a');
                    fprintf(fid, "'data_hist.pdf': file, the summary PDF of history-based data profiles for all tolerances.\n");
                    fprintf(fid, "'data_out.pdf': file, the summary PDF of output-based data profiles for all tolerances.\n");
                    fprintf(fid, "'perf_hist.pdf': file, the summary PDF of history-based performance profiles for all tolerances.\n");
                    fprintf(fid, "'perf_out.pdf': file, the summary PDF of output-based performance profiles for all tolerances.\n");
                    if n_solvers == 2
                        fprintf(fid, "'log-ratio_hist.pdf': file, the summary PDF of history-based log-ratio profiles for all tolerances.\n");
                        fprintf(fid, "'log-ratio_out.pdf': file, the summary PDF of output-based log-ratio profiles for all tolerances.\n");
                    end
                    fclose(fid);
                catch
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

    % Record the names of the problems all the solvers failed to meet the convergence test for every tolerance.
    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value) && (any(solvers_all_diverge_hist(:)) || any(solvers_all_diverge_out(:)))
        try
            fid = fopen(path_txt, 'a');
            fprintf(fid, "\n");
            fprintf(fid, "Problems that all the solvers failed to meet the convergence test for each tolerance and each run:\n");
            for i_tol = 1:profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)
                tolerance = tolerances(i_tol);
                tolerance_str = formatFloatScientificLatex(tolerance, 1);
                for i_run = 1:n_runs
                    if ~any(solvers_all_diverge_hist(:, i_run, i_tol)) && ~any(solvers_all_diverge_out(:, i_run, i_tol))
                        continue;
                    end
                    fprintf(fid, "History-based  tol = %-8s run = %-3d:\t\t", tolerance_str, i_run);
                    for i_problem = 1:n_problems
                        if solvers_all_diverge_hist(i_problem, i_run, i_tol)
                            fprintf(fid, "%-*s ", max_name_length, problem_names{i_problem});
                        end
                    end
                    fprintf(fid, "\n");
                    fprintf(fid, "Output-based   tol = %-8s run = %-3d:\t\t", tolerance_str, i_run);
                    for i_problem = 1:n_problems
                        if solvers_all_diverge_out(i_problem, i_run, i_tol)
                            fprintf(fid, "%-*s ", max_name_length, problem_names{i_problem});
                        end
                    end
                    fprintf(fid, "\n");
                end
            end
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to record the problems that all the solvers failed to meet the convergence test.\n");
            end
        end
    end

    if n_rows > 0 && ispc
        % Check the operating system. If it is Windows, we will adjust the position, papersize, paperposition of the
        % summary figure! (MATLAB in Windows will adjust the figure size automatically when it is too large.)
        fig_summary.Position = [defaultFigurePosition(1:2), profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * default_width, multiplier * n_rows * default_height];
        fig_summary.Units = 'centimeters';
        fig_summary.PaperUnits = 'centimeters';
        fig_summary.PaperSize = fig_summary.Position(3:4);
        fig_summary.PaperPosition = [0, 0, fig_summary.Position(3:4)];
        fig_summary.Children.Title.FontSize = min(fig_summary.Position(3) / 75 * 10, fig_summary.Position(3) / profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * 3 / 75 * 10);
        fig_summary.Children.Children(1).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
        fig_summary.Children.Children(2).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Save `curves` to path_log.
        save(fullfile(path_log, 'curves.mat'), 'curves');
        try
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'curves.mat': file, storing the curves of the profiles.\n");
            fclose(fid);
        catch
        end
    end

    % Compute the `solver_scores`.
    profile_scores = computeScores(curves, profile_options);
    score_fun = profile_options.(ProfileOptionKey.SCORE_FUN.value);
    solver_scores = score_fun(profile_scores);
    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        save(fullfile(path_log, 'profile_scores.mat'), 'profile_scores');
        try
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'profile_scores.mat': file, storing the scores of solvers on each profile.\n");
            fclose(fid);
        catch
        end
    end

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        % Print the scores of the solvers.
        fprintf('\n');
        fprintf('INFO: Scores of the solvers:\n');
        max_solver_name_length = max(cellfun(@length, solver_names));
        for i_solver = 1:n_solvers
            format_info_str = sprintf('INFO: %%-%ds:    %%.4f\n', max_solver_name_length);
            fprintf(format_info_str, solver_names{i_solver}, solver_scores(i_solver));
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Store the summary pdf. We will name the summary pdf as "summary_feature_name.pdf" and store it under
        % path_feature. We will also put a "summary.pdf" in the path_out directory, which will be a merged pdf of
        % all the "summary_feature_name.pdf" under path_out following the order of the feature_stamp.
        summary_name = ['summary_', feature_stamp, '_', time_stamp];
        if n_rows > 0
            if ispc
                print(fig_summary, fullfile(path_feature, [summary_name, '.pdf']), '-dpdf', '-vector');
            else
                exportgraphics(fig_summary, fullfile(path_feature, [summary_name, '.pdf']), 'ContentType', 'vector');
            end
            savefig(fig_summary, fullfile(path_figs, [summary_name, '.fig']));
        end

        try
            fid = fopen(path_readme_feature, 'a');
            fprintf(fid, "'%s': file, the summary PDF of all the profiles for the feature '%s'.\n", summary_name, feature_name);
            fclose(fid);
        catch
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
        if ~profile_options.(ProfileOptionKey.SILENT.value) && ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
            fprintf('\nINFO: Summary stored in: \n%s\n\n', path_out);
        end
    end

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('\nINFO: Finished the computation of the profiles under "%s" feature.\n', feature.name);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Save useful data to a structure and then to a mat file for future loading. %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        data = struct();
        data.feature_stamp = feature_stamp;
        data.solver_names = solver_names;
        data.n_evals = n_evals;
        data.problem_names = problem_names;
        data.problem_types = problem_types;
        data.problem_dims = problem_dims;
        data.problem_cons = problem_cons;
        data.computation_times = computation_times;
        data.solvers_successes = solvers_successes;
        data.merit_histories = merit_histories;
        data.merit_out = merit_out;
        data.merit_init = merit_init;
        try
            save(fullfile(path_log, 'data.mat'), '-struct', 'data');
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'data.mat': file, storing the data of the current experiment for future loading.\n");
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("INFO: Failed to save the data of the current experiment.\n");
            end
        end

        diary off;
    end

end

% Following one function is modified from the code provided by Benjamin Großmann (2024). Merge PDF-Documents
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

function integral = integrate(curve, profile_type, profile_options)
    % Compute the integral of the curve from different types of profiles.

    kernel = profile_options.(ProfileOptionKey.SCORE_WEIGHT_FUN.value);
    integral = 0;
    switch profile_type
        case {'perf', 'data'}
            % The curve is a right-continuous step function.
            integral = integral + sum(diff(curve(1, :)) .* curve(2, 1:end-1) .* kernel(curve(1, 1:end-1)));
        case 'log_ratio'
            % We do not modify the integral of log_ratio even a score_weight_fun is provided.
            if isempty(curve)
                integral = 0;
            else
                integral = integral + sum(abs(curve(2, :)));
            end
    end
end

function profile_scores = computeScores(curves, profile_options)
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
            profile_scores(i_solver, i_tol, 1, 1) = integrate(curve_hist_perf, 'perf', profile_options);
            profile_scores(i_solver, i_tol, 1, 2) = integrate(curve_hist_data, 'data', profile_options);
            profile_scores(i_solver, i_tol, 2, 1) = integrate(curve_out_perf, 'perf', profile_options);
            profile_scores(i_solver, i_tol, 2, 2) = integrate(curve_out_data, 'data', profile_options);
            if n_solvers == 2
                curve_hist_log_ratio = curves{i_tol}.hist.log_ratio{i_solver};
                curve_out_log_ratio = curves{i_tol}.out.log_ratio{i_solver};
                profile_scores(i_solver, i_tol, 1, 3) = integrate(curve_hist_log_ratio, 'log_ratio', profile_options);
                profile_scores(i_solver, i_tol, 2, 3) = integrate(curve_out_log_ratio, 'log_ratio', profile_options);
            end
        end
    end

    % Normalize the profile_scores by dividing the maximum profile_scores with the same tolerance and types.
    if profile_options.(ProfileOptionKey.NORMALIZED_SCORES.value)
        max_scores = max(profile_scores, [], 1);
        max_scores(max_scores == 0) = 1;
        profile_scores = profile_scores ./ max_scores;
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