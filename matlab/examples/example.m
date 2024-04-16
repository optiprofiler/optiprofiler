function example()

    clc;

    %%%%%%%%%%%%%%%%%%%%% custom_problem_names %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    folderPath = '~/Library/matlab_alg/opm/problems';  % The path to OPM dir.
    addpath(folderPath);
    excludeList = {'opm_eval_cpsf', 'shydc20ls', 'tcontact', 'deconvu', ...
        'eigenals', 'eigenbls', 'eigencls', 'hydc20ls', 'nondquar'};
    custom_problem_names = getMFileNames(folderPath, excludeList);
    % custom_problem_names = {'arwhead', 'arglinb', 'bard'};

    %%%%%%%%%%%%%%%%%%%%% cutest_problem_names %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cutest_struct.type = 'u';
    cutest_struct.maxdim = 200;
    cutest_struct.mindim = 6;
    cutest_struct.blacklist = {'ARGTRIGLS', 'BROWNAL', 'COATING', ...
        'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', 'DMN15103LS', ...
        'DMN15332LS', 'DMN15333LS', 'DMN37142LS', 'DMN37143LS', ...
        'ERRINRSM', 'HYDC20LS', 'LRA9A', 'LRCOVTYPE', 'LUKSAN12LS', ...
        'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', 'LUKSAN22LS', ...
        'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM'};
    cutest_problem_names = secup(cutest_struct);

    %%%%%%%%%%%%%%%%%%%%% remaining parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    solvers = {@fminsearch_test, @bds_test};
    labels = {'fminsearch', 'bds',};

    runBenchmark(solvers, labels, cutest_problem_names, @OPM_loader, custom_problem_names, 'all', 'max_tol_order', 10, 'summarize_log_ratio_profiles', true, 'benchmark_id', 'test_for_slides')

end

function x = fminsearch_test(fun, x0)
  
    x = fminsearch(fun, x0);
    
end

function x = bds_test(fun, x0)

    x = bds(fun, x0);

end