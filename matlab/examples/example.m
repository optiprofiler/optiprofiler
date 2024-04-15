function example()

    clc;
    % folderPath = '/Users/huangcunxin/Bureau/OPM/problems';
    folderPath = '~/Library/matlab_alg/opm/problems';  % The path to OPM dir.
    excludeList = {'opm_eval_cpsf'};
    custom_problem_names = getMFileNames(folderPath, excludeList);
    % custom_problem_names = {'arwhead', 'arglinb', 'bard'};

    cutest_problem_names = secup(struct('type', 'u', 'maxdim', 3));
    solvers = {@fminsearch_test, @bds_test};
    labels = {'fminsearch', 'bds',};

    runBenchmark(solvers, labels, cutest_problem_names, @OPM_loader, custom_problem_names, 'plain', 'n_runs', 1, 'max_tol_order', 10, 'summarize_log_ratio_profiles', true, 'benchmark_id', 'test_custom_problems')

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