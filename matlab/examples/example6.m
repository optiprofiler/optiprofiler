function example6()
%EXAMPLE6 corresponds to "Example 6: wrapping a solver with nonlinear constraints" in the
% "Usage for MATLAB" part of our official website (www.optprof.com).
%
% This example shows how to wrap a solver such as fmincon. OptiProfiler
% provides nonlinear inequalities and equalities as two callbacks, cub and
% ceq, while fmincon expects one callback with two outputs. The wrapper uses
% deal to perform this signature conversion.

    if exist('fmincon', 'file') ~= 2
        error('OptiProfiler:Example6:MissingFmincon', ...
            'This example requires fmincon from MATLAB Optimization Toolbox.');
    end

    % Print the information about this example.
    fprintf('\nThis example benchmarks fmincon through OptiProfiler wrappers.\n');
    fprintf('The wrappers combine OptiProfiler''s cub and ceq callbacks into one nonlcon callback.\n');
    pause(1.5);
    fprintf('\nStart Example 6...\n\n');

    % Select small constrained problems so that the example finishes quickly.
    options.ptype = 'n';
    options.problem_names = {'HS10', 'HS11', 'HS12'};
    options.plibs = {'s2mpj'};
    options.mindim = 2;
    options.maxdim = 5;
    options.max_eval_factor = 500;
    options.draw_hist_plots = 'none';
    options.n_jobs = 1;
    options.solver_names = {'fmincon short', 'fmincon long'};

    scores = benchmark({@fmincon_short, @fmincon_long}, options)
end

function x = fmincon_short(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)
%FMINCON_SHORT runs fmincon with a small function-evaluation budget.

    x = fmincon_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, 100);
end

function x = fmincon_long(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)
%FMINCON_LONG runs fmincon with a larger function-evaluation budget.

    x = fmincon_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, 200);
end

function x = fmincon_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_fun_evals)
%FMINCON_WRAPPER adapts OptiProfiler's nonlinear-constraint signature to fmincon.
%
% OptiProfiler provides nonlinear inequalities and equalities separately as
% cub(x) <= 0 and ceq(x) = 0. MATLAB solvers such as fmincon expect a single
% nonlcon callback returning both values: [c, ceq] = nonlcon(x). The call to
% deal below evaluates cub(x) and ceq(x), then returns them as those two
% outputs in the order fmincon expects.

    nonlcon = @(x) deal(cub(x), ceq(x));
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', max_fun_evals);
    x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
end
