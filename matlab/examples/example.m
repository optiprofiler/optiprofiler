function example()

    clear;
    clc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@bds_test, @fminsearch_test};
    % benchmark(solvers)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@bds_test, @fminsearch_test};
    % benchmark(solvers, 'noisy')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solvers = {@bds_test, @fminsearch_test};
    options.feature_names = {'plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted', };
    % options.noise_level = 1e-1;
    options.n_jobs = 1;
    options.cutest_problem_names = 'SISSER';
    % options.rotated = true;
    % options.condition_number = 'dimension_dependent';
    options.run_plain = false;
    % options.n_runs = 5;
    % options.max_tol_order = 4;
    % options.problem_type = 'u';
    % options.maxdim = 5;
    % options.benchmark_id = 'test-path';
    % options.summarize_log_ratio_profiles = true;
    options.summarize_maxcvhist = false;
    options.summarize_merithist = false;
    options.labels = {'bds', 'simplex'};
    benchmark(solvers, options)

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