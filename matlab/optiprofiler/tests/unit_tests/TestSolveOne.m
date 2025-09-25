classdef TestSolveOne < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            problem_name = 'ALLINITU';
            problem = s2mpj_load(problem_name);
            solvers = {@fminsearch_test1, @fminsearch_test2};
            feature = Feature('permuted', 'n_runs', 1);
            len_problem_names = length('ALLINITU');
            profile_options.solver_names = {'fminsearch', 'fminunc'};
            profile_options.solver_isrand = [false, false];
            profile_options.project_x0 = true;
            profile_options.max_eval_factor = 500;
            profile_options.silent = false;
            profile_options.seed = 1;
            profile_options.solver_verbose = 1;
            profile_options.score_only = false;
            profile_options = getDefaultProfileOptions(solvers, feature, profile_options);

            result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, true, '');
            testCase.verifyNotEmpty(result);
            testCase.verifyNotEmpty(result.fun_history);
            testCase.verifyNotEmpty(result.maxcv_history);
            testCase.verifyNotEmpty(result.fun_out);
            testCase.verifyNotEmpty(result.maxcv_out);
            testCase.verifyNotEmpty(result.fun_inits);
            testCase.verifyNotEmpty(result.maxcv_inits);
            testCase.verifyNotEmpty(result.n_eval);
            testCase.verifyNotEmpty(result.problem_name);
            testCase.verifyNotEmpty(result.problem_type);
            testCase.verifyNotEmpty(result.problem_dim);
            testCase.verifyNotEmpty(result.problem_mb);
            testCase.verifyNotEmpty(result.problem_mcon);
            testCase.verifyNotEmpty(result.computation_time);
            testCase.verifyNotEmpty(result.solvers_success);

            problem_name = 'BT1';
            problem = s2mpj_load(problem_name);
            solvers = {@fmincon_test1, @fmincon_test2};
            len_problem_names = length('BT1');
            profile_options.solver_names = {'sqp', 'interior-point'};
            result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, true, '');
            testCase.verifyNotEmpty(result);
            testCase.verifyNotEmpty(result.fun_history);
            testCase.verifyNotEmpty(result.maxcv_history);
            testCase.verifyNotEmpty(result.fun_out);
            testCase.verifyNotEmpty(result.maxcv_out);
            testCase.verifyNotEmpty(result.fun_inits);
            testCase.verifyNotEmpty(result.maxcv_inits);
            testCase.verifyNotEmpty(result.n_eval);
            testCase.verifyNotEmpty(result.problem_name);
            testCase.verifyNotEmpty(result.problem_type);
            testCase.verifyNotEmpty(result.problem_dim);
            testCase.verifyNotEmpty(result.problem_mb);
            testCase.verifyNotEmpty(result.problem_mcon);
            testCase.verifyNotEmpty(result.computation_time);
            testCase.verifyNotEmpty(result.solvers_success);
        end
    end
end

function x = fminsearch_test1(fun, x0)

x = fminsearch(fun, x0);
end

function x = fminsearch_test2(fun, x0)

options = optimset('MaxFunEvals', 1000);
x = fminunc(fun, x0, options);
end

function x = fmincon_test1(varargin)
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'sqp');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end

function x = fmincon_test2(varargin)

    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'interior-point');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end