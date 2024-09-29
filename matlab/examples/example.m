function example()

    clear;
    clc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solvers = {@bds_test, @fminsearch_test};
    benchmark(solvers)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@bds_test, @fminsearch_test};
    % benchmark(solvers, 'noisy')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@bds_test, @fminsearch_test};
    % feature_name = 'noisy';
    % problem = loader('TOINTGOR');
    % benchmark(solvers, feature_name, problem)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@bds_test, @fminsearch_test};
    % options.feature_name = 'plain';
    % options.benchmark_id = 'test';
    % options.run_plain = false;
    % options.problem_type = 'u';
    % options.maxdim = 3;
    % options.summarize_log_ratio_profiles = true;
    % options.labels = {'bds', 'simplex'};
    % benchmark(solvers, options)

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