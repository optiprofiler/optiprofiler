function example()

    addpath("src/")
    addpath(genpath("common/"))

    problem_names = secup(struct('type', 'u', 'maxdim', 3, 'mindim', 1));
    solvers = {@newuoa_test, @uobyqa_test};
    labels = {"newuoa", 'uobyqa'};

    runBenchmarks(solvers, labels, problem_names, 'all', 'benchmark_id', 'test_all_features', 'max_tol_order', 10)

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

function fminsearch_test(fun, x0)
  
    x = fminsearch(fun, x0);
    
end

function fminunc(fun, x0)
  
    x = fminunc(fun, x0);
    
end