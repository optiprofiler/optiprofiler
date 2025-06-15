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
%   SOLVER_SCORES = BENCHMARK(SOLVERS, FEATURE_NAME) creates profiles for the
%   given SOLVERS on the default unconstrained problem set with the specified
%   feature FEATURE_NAME.
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
%   Cautions:
%
%   1. Each solver in SOLVERS should accept corresponding signature(s)
%      depending on the test suite you choose:
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
%      solvers. For more information of performance and data profiles, see
%      [1]_, [2]_, [5]_. For that of log-ratio profiles, see [4]_, [6]_.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Options:
%
%   Options should be specified in a struct. The following are the available
%   fields of the struct:
%
%       1. Options for profiles and plots:
%
%       - n_jobs: the number of parallel jobs to run the test. Default is the
%         default number of workers in the default local cluster.
%       - seed: the seed of the random number generator. Default is 1.
%       - benchmark_id: the identifier of the test. It is used to create the
%         specific directory to store the results. Default is 'out' if the
%         option `load` is not provided, otherwise default is '.'.
%       - solver_names: the names of the solvers. Default is the names of the
%         function handles in `solvers`.
%       - solver_isrand: whether the solvers are randomized or not. Default is
%         a logical array of the same length as the number of solvers, where
%         the value is true if the solver is randomized, and false otherwise.
%         Note that if `n_runs` is not specified, we will set it 5 for the
%         randomized solvers.
%       - feature_stamp: the stamp of the feature with the given options. It is
%         used to create the specific directory to store the results. Default
%         depends on features.
%       - errorbar_type: the type of the uncertainty interval that can be
%         either 'minmax' or 'meanstd'. When `n_runs` is greater than 1, we run
%         several times of the experiments and get average curves and
%         uncertainty intervals. Default is 'minmax', meaning that we takes the
%         pointwise minimum and maximum of the curves.
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
%               ``(fun_value, maxcv_value, maxcv_init) -> merit_value``,
%         where `fun_value` is the objective function value, `maxcv_value` is
%         the maximum constraint violation, and `maxcv_init` is the maximum
%         constraint violation at the initial guess. The default merit function
%         varphi(x) is defined by the objective function f(x) and the maximum
%         constraint violation v(x) as
%           varphi(x) = f(x)                        if v(x) <= v1
%           varphi(x) = f(x) + 1e5 * (v(x) - v1)    if v1 < v(x) <= v2
%           varphi(x) = Inf                         if v(x) > v2
%         where v1 = max(1e-5, v0) and v2 = min(0.01, 1e-10 * max(1, v0)),
%         and v0 is the maximum constraint violation at the initial guess.
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
%       - solvers_to_load: the indices of the solvers to load when the `load`
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
%         '--', ':') and markers ('o', '+', '*', '.', 'x', 's', 'd', '^', 'v',
%         '>', '<', 'p', 'h'). Default line style order is {'-', '-.', '--',
%         ':'}. Note that if the number of solvers is greater than the number
%         of line styles, we will cycle through the styles.
%       - line_widths: the widths of the lines in the plots. It should be a
%         positive scalar or a vector. Default is 1.5. Note that if the number
%         of solvers is greater than the number of line widths, we will cycle
%         through the widths.
%       - bar_colors: two different colors for the bars of two solvers in the
%         log-ratio profiles. It can be a cell array of short names of colors
%         ('r', 'g', 'b', 'c', 'm', 'y', 'k') or a 2-by-3 matrix with each row
%         being a RGB triplet. Default is set to the first two colors in the
%         `line_colors` option.
%       - xlabel_performance_profile: the label of the x-axis of the
%         performance profiles. Default is 'Performance ratio'.
%         Note: the 'Interpreter' property is set to 'latex', so LaTeX
%         formatting is supported. The same applies to the options
%         `ylabel_performance_profile`, `xlabel_data_profile`,
%         `ylabel_data_profile`, `xlabel_log_ratio_profile`, and
%         `ylabel_log_ratio_profile`.
%       - ylabel_performance_profile: the label of the y-axis of the
%         performance profiles. Default is
%         'Performance profiles ($\\mathrm{tol} = %s$)', where `%s` will be
%         replaced by the current tolerance in LaTeX format. You can also use
%         `%s` in your custom label, and it will be replaced accordingly.
%       - xlabel_data_profile: the label of the x-axis of the data profiles.
%         Default is 'Number of simplex gradients'.
%       - ylabel_data_profile: the label of the y-axis of the data profiles.
%         Default is 'Data profiles ($\\mathrm{tol} = %s$)', where `%s` will be
%         replaced by the current tolerance in LaTeX format. You can also use
%         `%s` in your custom label, and it will be replaced accordingly.
%       - xlabel_log_ratio_profile: the label of the x-axis of the log-ratio
%         profiles. Default is 'Problem'.
%       - ylabel_log_ratio_profile: the label of the y-axis of the log-ratio
%         profiles. Default is 'Log-ratio profiles ($\\mathrm{tol} = %s$)',
%         where `%s` will be replaced by the current tolerance in LaTeX format.
%         You can also use `%s` in your custom label, and it will be replaced
%         accordingly.
%
%       2. Options for features:
%
%       - feature_name: the name of the feature. The available features are
%         'plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted',
%         'linearly_transformed', 'random_nan', 'unrelaxable_constraints',
%         'nonquantifiable_constraints', 'quantized', and 'custom'. Default is
%         'plain'.
%       - n_runs: the number of runs of the experiments with the given feature.
%         Default is 5 for stochastic features and 1 for deterministic
%         features.
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
%       - perturbed_trailing_digits: whether we will randomize the trailing
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
%       3. Options for problems:
%
%       Options in this part are used to select problems for benchmarking.
%       First select which problem libraries to use based on the `plibs`
%       option. Then select problems from these libraries according to the
%       given options (`problem_names`, `ptype`, `mindim`, `maxdim`, `minb`,
%       `maxb`, `minlcon`, `maxlcon`, `minnlcon`, `maxnlcon`, `mincon`,
%       `maxcon`, and `excludelist`).
%
%       Following is the list of available options:
%
%       - plibs: the problem libraries to be used. It should be a cell array of
%         strings or chars. The available choices are subfolder names in the
%         'problems' directory. There are three subfolders after installing the
%         package: 's2mpj', 'matcutest', and 'custom'. Default setting is
%         's2mpj'.
%       - ptype: the type of the problems to be selected. It should be a string
%         or char consisting of any combination of 'u' (unconstrained), 'b'
%         (bound constrained), 'l' (linearly constrained), and 'n' (nonlinearly
%         constrained), such as 'b', 'ul', 'ubn'. Default is 'u'.
%       - mindim: the minimum dimension of the problems to be selected. Default
%         is 1.
%       - maxdim: the maximum dimension of the problems to be selected. Default
%         is mindim + 1.
%       - minb: the minimum number of bound constraints of the problems to be
%         selected. Default is 0.
%       - maxb: the maximum number of bound constraints of the problems to be
%         selected. Default is minb + 10.
%       - minlcon: the minimum number of linear constraints of the problems to
%         be selected. Default is 0.
%       - maxlcon: the maximum number of linear constraints of the problems to
%         be selected. Default is minlcon + 10.
%       - minnlcon: the minimum number of nonlinear constraints of the problems
%         to be selected. Default is 0.
%       - maxnlcon: the maximum number of nonlinear constraints of the problems
%         to be selected. Default is minnlcon + 10.
%       - mincon: the minimum number of linear and nonlinear constraints of the
%         problems to be selected. Default is min(minlcon, minnlcon).
%       - maxcon: the maximum number of linear and nonlinear constraints of the
%         problems to be selected. Default is max(maxlcon, maxnlcon).
%       - excludelist: the list of problems to be excluded. Default is not to
%         exclude any problem.
%       - problem_names: the names of the problems to be selected. It should
%         be a cell array of strings or chars. Default is not to select any
%         problem by name but by the options above.
%
%       You may also pass an instance of the class Problem by the option
%
%       - problem: a problem to be benchmarked. It should be an instance of the
%         class Problem. If it is provided, we will only run the test on this
%         problem with the given feature and draw the history plots. Default is
%         not to set any problem.
%
%   Several points to note:
%
%   1. The information about two problem libraries is available in the
%      following links:
%           S2MPJ (see [3]_) <https://github.com/GrattonToint/S2MPJ>
%           MatCUTEst <https://github.com/matcutest>
%
%   2. If you want to use your own problem library, please check the README.txt
%      in the directory 'problems/' or the guidance in our website
%      <https://optprof.com> for more details.
%
%   3. The problem library MatCUTEst is only available when the OS is Linux.
%
%   4. If the `load` option is provided, we will use the provided options to
%      select data from the specified experiment for plotting the profiles.
%      Available options are:
%      - Options for profiles and plots: `benchmark_id`, `solver_names`,
%        `feature_stamp`, `errorbar_type`, `savepath`, `max_tol_order`,
%        `merit_fun`, `run_plain`, `score_only`,
%        `summarize_performance_profiles`, `summarize_data_profiles`,
%        `summarize_log_ratio_profiles`, `summarize_output_based_profiles`,
%        `silent`, `semilogx`, `normalized_scores`, `score_weight_fun`,
%        `score_fun`, `solvers_to_load`, `line_colors`, `line_styles`,
%        `line_widths`, `bar_colors`.
%      - Options for features: none.
%      - Options for problems: `plibs`, `ptype`, `mindim`, `maxdim`, `minb`,
%        `maxb`, `minlcon`, `maxlcon`, `minnlcon`, `maxnlcon`, `mincon`,
%        `maxcon`, `excludelist`.
%
%   5. More information about OptiProfiler can be found on our website:
%
%                           https://www.optprof.com
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
        error("MATLAB:benchmark:solverMustBeProvided", "At least a cell of function handles (callable solvers) or a struct of options must be provided to `benchmark`.");
    elseif nargin == 1
        if isstruct(varargin{1})
            % When input contains one argument and the first argument is a struct, we assume the
            % user chooses benchmark(options).
            options = varargin{1};
            options_user = options;
            solvers = {};
            if isfield(options, 'feature_name')
                feature_name = options.feature_name;
                options = rmfield(options, 'feature_name');
            else
                feature_name = 'plain';
            end
            if isfield(options, 'problem')
                options = rmfield(options, 'problem');
            end
            if ~isfield(options, ProfileOptionKey.LOAD.value) || isempty(options.(ProfileOptionKey.LOAD.value))
                error("MATLAB:benchmark:LoadFieldNotProvided", "The field `load` must be provided when the first argument for `benchmark` is a struct.");
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
            options_user = options;
            if isfield(options, 'feature_name')
                feature_name = options.feature_name;
                options = rmfield(options, 'feature_name');
            else
                feature_name = 'plain';
            end
            if isfield(options, 'problem')
                problem = options.problem;
                options = rmfield(options, 'problem');
            end
        else
            error("MATLAB:benchmark:SecondArgumentWrongType", "The second argument for `benchmark` must be a feature name or a struct of options.");
        end
    else
        error("MATLAB:benchmark:TooMuchInput", "Invalid number of arguments. The function `benchmark` at most takes two arguments.");
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
        error("MATLAB:benchmark:feature_nameNotcharstr", "`feature_name` provided for `benchmark` must be a char or string.");
    end
    feature_name = char(lower(feature_name));
    valid_feature_names = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
    if ~ismember(feature_name, valid_feature_names)
        error("MATLAB:benchmark:feature_nameNotValid", "`feature_name` provided for `benchmark` must be one of the valid feature names: %s.", strjoin(valid_feature_names, ', '));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Parse options for feature, problem, and profile %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    feature_options = struct();
    problem_options = struct();
    profile_options = struct();

    fieldNames = fieldnames(options);
    for i_field = 1:numel(fieldNames)
        key = fieldNames{i_field};
        value = options.(key);

        % We can also use: validFeatureOptionKeys = {enumeration('FeatureOptionKey').value};  
        % However, it only works for MATLAB R2021b or later.
        validFeatureOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('FeatureOptionKey')), 'UniformOutput', false);
        validProblemOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('ProblemOptionKey')), 'UniformOutput', false);
        validProfileOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('ProfileOptionKey')), 'UniformOutput', false);

        if ismember(key, validFeatureOptionKeys)
            feature_options.(key) = value;
        elseif ismember(key, validProblemOptionKeys)
            problem_options.(key) = value;
        elseif ismember(key, validProfileOptionKeys)
            profile_options.(key) = value;
        else
            error("MATLAB:benchmark:UnknownOptions", "Unknown option for `benchmark`: %s", key);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check validity of options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note: the validity of the feature options has been checked in the Feature constructor, so we
    % do not need to check it here.
    problem_options = checkValidityProblemOptions(problem_options, profile_options);
    profile_options = checkValidityProfileOptions(solvers, profile_options);

    % Whether load the existing results.
    is_load = isfield(profile_options, ProfileOptionKey.LOAD.value) && ~isempty(profile_options.(ProfileOptionKey.LOAD.value));

    % If `n_runs` is not specified, we set it to 5 if at least one solver is randomized.
    any_solver_isrand = isfield(profile_options, ProfileOptionKey.SOLVER_ISRAND.value) && any(profile_options.(ProfileOptionKey.SOLVER_ISRAND.value));
    if ~isfield(feature_options, FeatureOptionKey.N_RUNS.value) && any_solver_isrand && ~is_load
        if ~isfield(profile_options, ProfileOptionKey.SILENT.value) || ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf("\nINFO: We set `n_runs` to 5 since it is not specified and at least one solver is randomized.\n");
        end
        feature_options.(FeatureOptionKey.N_RUNS.value) = 5;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Process the 'load' option if it is provided. %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [results_plibs, profile_options] = loadResults(problem_options, profile_options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get default options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Build feature.
    feature = Feature(feature_name, feature_options);
    feature_options = feature.options;

    % Set default values for the unspecified options.
    problem_options = getDefaultProblemOptions(problem_options);
    profile_options = getDefaultProfileOptions(solvers, feature, profile_options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize output variables. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist('results_plibs', 'var') && ~isempty(results_plibs)
        n_solvers = size(results_plibs{1}.fun_histories, 2);
    elseif isfield(profile_options, ProfileOptionKey.SOLVER_NAMES.value)
        n_solvers = numel(profile_options.(ProfileOptionKey.SOLVER_NAMES.value));
    else
        n_solvers = numel(solvers);
    end
    solver_scores = zeros(n_solvers, 1);
    profile_scores = [];
    curves = [];
    solver_names = profile_options.(ProfileOptionKey.SOLVER_NAMES.value);
    % Remove backslash from the solver names.
    solver_names = cellfun(@(s) strrep(s, '\_', '_'), solver_names, 'UniformOutput', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Set the default options for plotting. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontName', 'Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Define the directory to store results. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define the directory to store or load the results.
    path_out = fullfile(profile_options.(ProfileOptionKey.SAVEPATH.value), profile_options.(ProfileOptionKey.BENCHMARK_ID.value));

    % Record the time stamp of the current experiment.
    time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
    time_stamp = char(time);

    % Set the feature stamp.
    feature_stamp = profile_options.(ProfileOptionKey.FEATURE_STAMP.value);

    % Create the stamp for the current experiment.
    stamp = createStamp(solver_names, problem_options, feature_stamp, time_stamp);

    path_stamp = fullfile(path_out, stamp);
    path_log = fullfile(path_stamp, 'test_log');
    path_report = fullfile(path_log, 'report.txt');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Create the directory to store the results. %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        if ~exist(path_stamp, 'dir')
            mkdir(path_stamp);
        else
            rmdir(path_stamp, 's');
            mkdir(path_stamp);
        end
    end

    if profile_options.(ProfileOptionKey.SCORE_ONLY.value) || is_load
        % If 'load' is not empty, we will directly load the results and do not need to compute
        % again. In this case, we do not need to create the directory to store the hist plots.
        path_hist_plots = '';
    else
        path_hist_plots = fullfile(path_stamp, 'history_plots');
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
            fprintf(fid, "# Content of this folder\n\n");
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("\nINFO: Failed to create the README.txt file for folder '%s'.\n", path_log);
            end
        end
    else
        path_readme_log = '';
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
                fprintf("\nINFO: Failed to create the time_stamp file for folder '%s'.\n", path_log);
            end
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value) && ~exist('problem', 'var')
        path_figs = fullfile(path_log, 'profile_figs');
        mkdir(path_figs);
        try
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'profile_figs': folder, containing all the FIG files of the profiles.\n");
            fclose(fid);
        catch
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        try
            if exist('options_user', 'var')
                save(fullfile(path_log, 'options_user.mat'), 'options_user');
                try
                    fid = fopen(path_readme_log, 'a');
                    fprintf(fid, "'options_user.mat': file, storing the options provided by the user for the current experiment.\n");
                    fclose(fid);
                catch
                end
                options_refined = profile_options;
                feature_options_keys = fieldnames(feature_options);
                for i_option = 1:numel(feature_options_keys)
                    key = feature_options_keys{i_option};
                    options_refined.(key) = feature_options.(key);
                end
                problem_options_keys = fieldnames(problem_options);
                for i_option = 1:numel(problem_options_keys)
                    key = problem_options_keys{i_option};
                    options_refined.(key) = problem_options.(key);
                end
                save(fullfile(path_log, 'options_refined.mat'), 'options_refined');
                try
                    fid = fopen(path_readme_log, 'a');
                    fprintf(fid, "'options_refined.mat': file, storing the options refined by OptiProfiler for the current experiment.\n");
                    fclose(fid);
                catch
                end
            end
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("\nINFO: Failed to save the options of the current experiment.\n");
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
        % Create a README.txt file to explain the content of the folder `path_stamp`.
        path_readme_feature = fullfile(path_stamp, 'README.txt');
        try
            fid = fopen(path_readme_feature, 'w');
            fprintf(fid, "# Content of this folder\n\n");
            if ~isempty(path_hist_plots)
                fprintf(fid, "'history_plots': folder, containing all the history plots for each problem.\n");
                fprintf(fid, "'history_plots_summary.pdf': file, the summary PDF of history plots for all problems.\n");
            end
            fprintf(fid, "'test_log': folder, containing log files and other useful experimental data.\n");
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("\nINFO: Failed to create the README.txt file for folder '%s'.\n", path_stamp);
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
                    fprintf("\nINFO: The script or function that calls `benchmark` function is copied to: %s.\n", path_log);
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
                fprintf("\nINFO: Failed to copy the script or function that calls `benchmark` function to the log directory.\n");
            end
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value) && ~exist('problem', 'var')
        % Create the directories for the performance profiles, data profiles, and log-ratio profiles.
        path_perf_hist = fullfile(path_stamp, 'detailed_profiles', 'perf_history-based');
        path_data_hist = fullfile(path_stamp, 'detailed_profiles', 'data_history-based');
        path_log_ratio_hist = fullfile(path_stamp, 'detailed_profiles', 'log-ratio_history-based');
        path_perf_out = fullfile(path_stamp, 'detailed_profiles', 'perf_output-based');
        path_data_out = fullfile(path_stamp, 'detailed_profiles', 'data_output-based');
        path_log_ratio_out = fullfile(path_stamp, 'detailed_profiles', 'log-ratio_output-based');
        
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

        pdf_perf_hist_summary = fullfile(path_stamp, 'perf_hist.pdf');
        pdf_perf_out_summary = fullfile(path_stamp, 'perf_out.pdf');
        pdf_data_hist_summary = fullfile(path_stamp, 'data_hist.pdf');
        pdf_data_out_summary = fullfile(path_stamp, 'data_out.pdf');
        pdf_log_ratio_hist_summary = fullfile(path_stamp, 'log-ratio_hist.pdf');
        pdf_log_ratio_out_summary = fullfile(path_stamp, 'log-ratio_out.pdf');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Solve all the problems when 'load' option is not provided. %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % If a specific problem is provided to `problem_options`, we only solve this problem and generate
    % the history plots for it.
    if exist('problem', 'var')
        result = solveOneProblem(solvers, problem, feature, problem.name, length(problem.name), profile_options, true, path_hist_plots);
        if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
            % We move the history plots to the feature directory.
            try
                movefile(fullfile(path_hist_plots, '*'), path_stamp);
                rmdir(path_hist_plots, 's');
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fprintf('\nINFO: Detailed results stored in: \n%s\n', path_stamp);
                end
            catch
            end
        end
        % Compute merit values.
        merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
        if ~isempty(result)
            try
                merit_history = meritFunCompute(merit_fun, result.fun_history, result.maxcv_history, result.maxcv_init);
                merit_init = meritFunCompute(merit_fun, result.fun_init, result.maxcv_init, result.maxcv_init);
            catch
                error("MATLAB:benchmark:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
            end
            % Find the least merit value for each problem.
            merit_min = min(merit_history, [], 'all');
            merit_min = min(merit_min, merit_init, 'omitnan');
            % Since we will not compute the profiles, we set `solver_scores` to be the relative decreases
            % in the objective function value.
            solver_merit_mins = squeeze(min(min(merit_history, [], 3, 'omitnan'), [], 2, 'omitnan'));
            solver_scores = (merit_init - solver_merit_mins) ./ max(merit_init - merit_min, eps);
        else
            solver_scores = zeros(n_solvers, 1);
        end

        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('\nINFO: Scores of the solvers:\n');
            max_solver_name_length = max(cellfun(@length, solver_names));
            for i_solver = 1:n_solvers
                format_info_str = sprintf('INFO: %%-%ds:    %%.4f\n', max_solver_name_length);
                fprintf(format_info_str, solver_names{i_solver}, solver_scores(i_solver));
            end
        end
        if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
            diary off;
        end
        return;
    end

    % If `load` is not specified, we solve all the problems.
    if ~is_load
        % Print the information about the experiment.
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('\nINFO: Start testing with the following options:\n');
            fprintf('INFO: - Solvers: %s\n', strjoin(solver_names, ', '));
            fprintf('INFO: - Problem libraries: %s\n', strjoin(problem_options.(ProblemOptionKey.PLIBS.value), ', '));
            fprintf('INFO: - Problem types: %s\n', problem_options.(ProblemOptionKey.PTYPE.value));
            fprintf('INFO: - Problem dimension range: [%d, %d]\n', problem_options.(ProblemOptionKey.MINDIM.value), problem_options.(ProblemOptionKey.MAXDIM.value));
            if any(ismember(problem_options.(ProblemOptionKey.PTYPE.value), 'bln'))
            fprintf('INFO: - Problem mb range: [%d, %d]\n', problem_options.(ProblemOptionKey.MINB.value), problem_options.(ProblemOptionKey.MAXB.value));
            end
            if any(ismember(problem_options.(ProblemOptionKey.PTYPE.value), 'ln'))
            fprintf('INFO: - Problem mlcon range: [%d, %d]\n', problem_options.(ProblemOptionKey.MINLCON.value), problem_options.(ProblemOptionKey.MAXLCON.value));
            end
            if any(ismember(problem_options.(ProblemOptionKey.PTYPE.value), 'n'))
            fprintf('INFO: - Problem mnlcon range: [%d, %d]\n', problem_options.(ProblemOptionKey.MINNLCON.value), problem_options.(ProblemOptionKey.MAXNLCON.value));
            end
            if any(ismember(problem_options.(ProblemOptionKey.PTYPE.value), 'ln'))
            fprintf('INFO: - Problem mcon range: [%d, %d]\n', problem_options.(ProblemOptionKey.MINCON.value), problem_options.(ProblemOptionKey.MAXCON.value));
            end
            if ~isempty(problem_options.(ProblemOptionKey.PROBLEM_NAMES.value))
            fprintf('INFO: - Number of user-provided problem names: %d\n', numel(problem_options.(ProblemOptionKey.PROBLEM_NAMES.value)));
            end
            if ~isempty(problem_options.(ProblemOptionKey.EXCLUDELIST.value))
            fprintf('INFO: - Number of user-excluded problem names: %d\n', numel(problem_options.(ProblemOptionKey.EXCLUDELIST.value)));
            end
            fprintf('INFO: - Feature stamp: %s\n', feature_stamp);
        end
        % We will solve all the problems from all the problem libraries that user specified in the `problem_options`.
        plibs = problem_options.(ProblemOptionKey.PLIBS.value);
        results_plibs = cell(1, numel(plibs));
        for i_plib = 1:numel(plibs)
            plib = plibs{i_plib};
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf('\nINFO: Start testing problems from the problem library "%s" with "%s" feature.\n', plib, feature.name);
                if strcmp(plib, 's2mpj')
                    fprintf("\nINFO: More information about the S2MPJ problem library can be found at: https://github.com/GrattonToint/S2MPJ\n");
                elseif strcmp(plib, 'matcutest')
                    fprintf("\nINFO: More information about the MatCUTEst problem library can be found at: https://github.com/matcutest\n");
                end
            end

            % Create directory to store the history plots for each problem library.
            if isempty(path_hist_plots)
                path_hist_plots_lib = '';
            else
                path_hist_plots_lib = fullfile(path_hist_plots, plib);
                if ~exist(path_hist_plots_lib, 'dir')
                    mkdir(path_hist_plots_lib);
                end
            end

            % Add the path of the problem library to the MATLAB path.          
            mydir = fileparts(mfilename('fullpath'));
            plib_path = fullfile(mydir, '../../../problems', plib);
            addpath(plib_path);

            % Solve all the problems from the current problem library with the specified options and
            % get the computation results.
            results_plib = solveAllProblems(solvers, plib, feature, problem_options, profile_options, true, path_hist_plots_lib);

            % If there are no problems selected or solved, skip the rest of the code, print a message,
            % and continue to the next library.
            if isempty(fieldnames(results_plib))
                results_plibs{i_plib} = {};
                continue;
            end

            % Compute merit values.
            merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
            try
                merit_histories = meritFunCompute(merit_fun, results_plib.fun_histories, results_plib.maxcv_histories, results_plib.maxcv_inits);
                merit_outs = meritFunCompute(merit_fun, results_plib.fun_outs, results_plib.maxcv_outs, results_plib.maxcv_inits);
                merit_inits = meritFunCompute(merit_fun, results_plib.fun_inits, results_plib.maxcv_inits, results_plib.maxcv_inits);
            catch
                error("MATLAB:benchmark:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
            end
            results_plib.merit_histories = merit_histories;
            results_plib.merit_outs = merit_outs;
            results_plib.merit_inits = merit_inits;

            % Run the 'plain' feature if run_plain is true.
            if profile_options.(ProfileOptionKey.RUN_PLAIN.value)
                feature_plain = Feature(FeatureName.PLAIN.value);
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fprintf('\nINFO: Start testing problems from the problem library "%s" with "plain" feature.\n', plib);
                end
                results_plib_plain = solveAllProblems(solvers, plib, feature_plain, problem_options, profile_options, false, {});
                try
                    results_plib_plain.merit_histories = meritFunCompute(merit_fun, results_plib_plain.fun_histories, results_plib_plain.maxcv_histories, results_plib_plain.maxcv_inits);
                    results_plib_plain.merit_outs = meritFunCompute(merit_fun, results_plib_plain.fun_outs, results_plib_plain.maxcv_outs, results_plib_plain.maxcv_inits);
                    results_plib_plain.merit_inits= meritFunCompute(merit_fun, results_plib_plain.fun_inits, results_plib_plain.maxcv_inits, results_plib_plain.maxcv_inits);
                catch
                    error("MATLAB:benchmark:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
                end

                % Store data of the 'plain' feature for later calculating merit_mins.
                results_plib.results_plib_plain = results_plib_plain;
            end
            
            results_plibs{i_plib} = results_plib;

            % Merge the history plots for each problem library to a single pdf file.
            if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value) && any(results_plib.solvers_successes(:))
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fprintf('\nINFO: Merging all the history plots of problems from the "%s" library to a single PDF file.\n', plib);
                end
                try
                    mergePdfs(path_hist_plots_lib, [plib '_history_plots_summary.pdf'], path_hist_plots);
                catch
                    if ~profile_options.(ProfileOptionKey.SILENT.value)
                        fprintf('\nINFO: Failed to merge the history plots to a single PDF file.\n');
                    end
                end
            end
        end
        % Remove the empty elements from the results_plibs cell array.
        results_plibs = results_plibs(~cellfun('isempty', results_plibs));
        % If results_plibs is empty, we will directly return.
        if isempty(results_plibs)
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("\nINFO: No problems are selected or solved from any problem library.\n");
            end
            return;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Store the data for loading. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        try
            save(fullfile(path_log, 'data_for_loading.mat'), 'results_plibs');
            fid = fopen(path_readme_log, 'a');
            fprintf(fid, "'data_for_loading.mat': file, storing the data of the current experiment for future loading.\n");
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("\nINFO: Failed to save the data of the current experiment.\n");
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write the report file. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    writeReport(profile_options, results_plibs, path_report, path_readme_log);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Start the computation of all profiles. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process the results from all the problem libraries.
    [merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged] = processResults(results_plibs, profile_options);

    [n_problems, n_solvers, n_runs, ~] = size(merit_histories_merged);

    if ~profile_options.(ProfileOptionKey.SILENT.value) 
        fprintf('\nINFO: Start computing profiles.\n');
    end

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
        T_stamp = strrep(stamp, '_', '\_');
        T_title = ['Profiles for ``', T_stamp, '"'];
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

        work_hist = NaN(n_problems, n_solvers, n_runs);
        work_out = NaN(n_problems, n_solvers, n_runs);
        for i_problem = 1:n_problems
            for i_solver = 1:n_solvers
                for i_run = 1:n_runs
                    if isfinite(merit_mins_merged(i_problem))
                        threshold = max(tolerance * merit_inits_merged(i_problem) + (1 - tolerance) * merit_mins_merged(i_problem), merit_mins_merged(i_problem));
                    else
                        threshold = -Inf;
                    end
                    if min(merit_histories_merged(i_problem, i_solver, i_run, :), [], 'omitnan') <= threshold
                        work_hist(i_problem, i_solver, i_run) = find(merit_histories_merged(i_problem, i_solver, i_run, :) <= threshold, 1, 'first');
                    end
                    if merit_outs_merged(i_problem, i_solver, i_run) <= threshold
                        work_out(i_problem, i_solver, i_run) = n_evals_merged(i_problem, i_solver, i_run);
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

        [fig_perf_hist, fig_data_hist, fig_log_ratio_hist, curves{i_tol}.hist] = drawProfiles(work_hist, problem_dims_merged, solver_names, tolerance_latex, cell_axs_summary_hist, true, is_perf, is_data, is_log_ratio, profile_options, curves{i_tol}.hist);
        [fig_perf_out, fig_data_out, fig_log_ratio_out, curves{i_tol}.out] = drawProfiles(work_out, problem_dims_merged, solver_names, tolerance_latex, cell_axs_summary_out, is_output_based, is_perf, is_data, is_log_ratio, profile_options, curves{i_tol}.out);

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
    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        try
            max_name_length = max(cellfun(@length, problem_names_merged));
            fid = fopen(path_report, 'a');
            fprintf(fid, "\n");
            fprintf(fid, "## Problems among all the libraries that all the solvers failed to meet the convergence test for each tolerance and each run\n");
            if (any(solvers_all_diverge_hist(:)) || any(solvers_all_diverge_out(:)))
                for i_tol = 1:profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)
                    tolerance = tolerances(i_tol);
                    tolerance_str = formatFloatScientificLatex(tolerance, 1);
                    for i_run = 1:n_runs
                        if any(solvers_all_diverge_hist(:, i_run, i_tol))
                            fprintf(fid, "\n");
                            fprintf(fid, "History-based  tol = %-8s run = %-3d:\t\t", tolerance_str, i_run);
                            for i_problem = 1:n_problems
                                if solvers_all_diverge_hist(i_problem, i_run, i_tol)
                                    fprintf(fid, "%-*s ", max_name_length, problem_names_merged{i_problem});
                                end
                            end
                        end
                        if any(solvers_all_diverge_out(:, i_run, i_tol))
                            fprintf(fid, "\n");
                            fprintf(fid, "Output-based   tol = %-8s run = %-3d:\t\t", tolerance_str, i_run);
                            for i_problem = 1:n_problems
                                if solvers_all_diverge_out(i_problem, i_run, i_tol)
                                    fprintf(fid, "%-*s ", max_name_length, problem_names_merged{i_problem});
                                end
                            end
                        end
                    end
                end
                fprintf(fid, "\n");
            else
                fprintf(fid, "\n");
                fprintf(fid, "This part is empty.\n");
            end
            fclose(fid);
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf("\nINFO: Failed to record the problems that all the solvers failed to meet the convergence test.\n");
            end
        end
    end

    if n_rows > 0 && ispc
        % Check the operating system. If it is Windows, we will adjust the position, papersize, paperposition
        % of the summary figure! (MATLAB in Windows will adjust the figure size automatically when it is too large.)
        fig_summary.Position = [defaultFigurePosition(1:2), profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * default_width, multiplier * n_rows * default_height];
        fig_summary.Units = 'centimeters';
        fig_summary.PaperUnits = 'centimeters';
        fig_summary.PaperSize = fig_summary.Position(3:4);
        fig_summary.PaperPosition = [0, 0, fig_summary.Position(3:4)];
        fig_summary.Children.Title.FontSize = min(fig_summary.Position(3) / 75 * 10, fig_summary.Position(3) / profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) * 3 / 75 * 10);
        if profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value)
            fig_summary.Children.Children(1).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
            fig_summary.Children.Children(2).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
        else
            fig_summary.Children.Children(1).YLabel.FontSize = min(fig_summary.Position(4) / 30 * 4.5, fig_summary.Position(4) / n_rows * 2 / 30 * 4.5);
        end
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

            fid = fopen(path_report, 'a');
            fprintf(fid, "\n");
            fprintf(fid, "## Scores of the solvers\n\n");
            max_solver_name_length = max(cellfun(@length, solver_names));
            for i_solver = 1:n_solvers
                fprintf(fid, "%-*s:    %.4f\n", max_solver_name_length, solver_names{i_solver}, solver_scores(i_solver));
            end
        catch
        end
    end

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        % Print the scores of the solvers.
        fprintf('\nINFO: Scores of the solvers:\n');
        max_solver_name_length = max(cellfun(@length, solver_names));
        for i_solver = 1:n_solvers
            format_info_str = sprintf('INFO: %%-%ds:    %%.4f\n', max_solver_name_length);
            fprintf(format_info_str, solver_names{i_solver}, solver_scores(i_solver));
        end
    end

    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        % Store the summary pdf. We will name the summary pdf as "summary_stamp.pdf" and store
        % it under path_stamp. We will also put a "summary.pdf" in the path_out directory, which will
        % be a merged pdf of latest 10 "summary_stamp.pdf" under path_out following descending order
        % of their time stamps.
        summary_name = ['summary_', stamp];
        if n_rows > 0
            if ispc
                print(fig_summary, fullfile(path_stamp, [summary_name, '.pdf']), '-dpdf', '-vector');
            else
                exportgraphics(fig_summary, fullfile(path_stamp, [summary_name, '.pdf']), 'ContentType', 'vector');
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
        % Sort the summary PDF files by their time stamps (last 15 characters of the file names).
        summary_time_stamps = arrayfun(@(f) datetime(f.name(end-18:end-4), 'InputFormat', 'yyyyMMdd_HHmmss'), summary_files);
        [~, idx] = sort(summary_time_stamps, 'descend');
        summary_files = summary_files(idx);
        % Keep only the latest 10 summary PDF files.
        summary_files = summary_files(1:min(10, numel(summary_files)));
        try
            % Merge all the summary PDF files to a single PDF file.
            delete(fullfile(path_out, 'summary.pdf'));
            for i_file = 1:numel(summary_files)
                copyfile(fullfile(summary_files(i_file).folder, summary_files(i_file).name), path_out);
                mergePdfs(path_out, 'summary.pdf', path_out);
                delete(fullfile(path_out, summary_files(i_file).name));
            end
        catch
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                fprintf('\nINFO: Failed to merge the summary PDF files.\n');
            end
        end
    end
    warning('on');

    % Close the figures.
    if n_rows > 0
        close(fig_summary);
    end

    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('\nINFO: Finished computing profiles with "%s" feature.\n', feature.name);
        fprintf("\nINFO: Report of the experiment ('report.txt') is stored in: \n%s\n\n", path_log);
        if ~isempty(path_hist_plots)
            fprintf('\nINFO: PDF of merged history plots is stored in: \n%s\n\n', path_hist_plots);
        end
        fprintf('\nINFO: Summary PDF of profiles is stored in: \n%s\n\n', path_out);
        fprintf('\nINFO: Single profiles are stored in: \n%s\n\n', fullfile(path_stamp, 'detailed_profiles'));
    end

    % Close the diary file.
    if ~profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        diary off;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunctions for the main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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