function example()

    clc;
    addpath(genpath('../optiprofiler/common'))

    problem_names = secup(struct('type', 'u', 'maxdim', 3));
    solvers = {@fminsearch_test, @fminunc_test, @fmincon_test};
    labels = {'fminsearch', 'fminunc', 'fmincon'};

    runBenchmark(solvers, labels, problem_names, 'plain', 'max_tol_order', 10, 'summarize_log_ratio_profiles', true, 'benchmark_id', 'test', 'n_jobs', 1)

end

function x = newuoa_test(fun, x0)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.solver = 'newuoa';

    x = pdfo(problem);

end

function x = uobyqa_test(fun, x0)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.solver = 'uobyqa';

    x = pdfo(problem);

end

function x = fminsearch_test(fun, x0)
  
    x = fminsearch(fun, x0);
    
end

function x = fminunc_test(fun, x0)
  
    x = fminunc(fun, x0);
    
end

function x = bds_test(fun, x0)

    x = bds(fun, x0);

end

function x = fmincon_test(fun, x0)

    x = fmincon(fun, x0);
end