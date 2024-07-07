function example()

    clc;

    %%%%%%%%%%%%%%%%%%%%% custom_problem_names %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % folderPath = '~/Library/matlab_alg/opm/problems';  % The path to OPM dir.
    % folderPath = '~/Bureau/OPM/problems';
    % addpath(folderPath);
    % excludeList = {'opm_eval_cpsf', 'shydc20ls', 'tcontact', 'deconvu', ...
    %     'eigenals', 'eigenbls', 'eigencls', 'hydc20ls', 'nondquar'};
    % custom_problem_names = getMFileNames(folderPath, excludeList);
    % custom_problem_names = {'arwhead', 'arglinb', 'bard'};

    %%%%%%%%%%%%%%%%%%%%% cutest_problem_names %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cutest_struct.type = 'u';
    % cutest_struct.maxdim = 200;
    % cutest_struct.mindim = 6;
    % cutest_struct.blacklist = {'ARGTRIGLS', 'BROWNAL', 'COATING', ...
    %     'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', 'DMN15103LS', ...
    %     'DMN15332LS', 'DMN15333LS', 'DMN37142LS', 'DMN37143LS', ...
    %     'ERRINRSM', 'HYDC20LS', 'LRA9A', 'LRCOVTYPE', 'LUKSAN12LS', ...
    %     'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', 'LUKSAN22LS', ...
    %     'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM'};
    % cutest_problem_names = {};

    %%%%%%%%%%%%%%%%%%%%% remaining parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    solvers = {@fminsearch_test, @bds_test};

    % benchmark(solvers, labels, {}, @OPM_loader, custom_problem_names, 'rotated', 'max_tol_order', 10, 'summarize_log_ratio_profiles', true, 'benchmark_id', 'test')
    % options.feature_names = 'plain';
    % options.max_tol_order = 10;
    % options.summarize_log_ratio_profiles = true;
    % options.n_jobs = 1;
    % options.labels = {"fminsearch", "BDS"};
    benchmark(solvers)

    % % % % % solvers = {@a, @a};
    % % % % % options.algorithm_options = {option1, option2};
    % % % % % benchmark(solvers, options);
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