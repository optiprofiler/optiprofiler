function test_u_plain()
    % Test `fminsearch` and `fminunc` with the "plain" feature.

    % Go to the directory of this repository.
    cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..'));

    solvers = {@fminsearch, @fminunc};
    options.solver_names = {'simplex', 'bfgs'};
    options.feature_name = 'plain';
    options.p_type = 'u';
    options.mindim = 1;
    options.maxdim = 2;
    options.benchmark_id = 'test_u_plain';

    x = benchmark(solvers, options)
end