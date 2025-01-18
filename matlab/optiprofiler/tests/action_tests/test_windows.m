function test_windows()

    % Go to the directory of this repository.
    cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..'));

    solvers = {@fminsearch, @fminunc};
    options.solver_names = {'simplex', 'bfgs'};
    options.p_type = 'u';
    options.mindim = 6;
    options.maxdim = 6;

    benchmark(solvers, options);
end