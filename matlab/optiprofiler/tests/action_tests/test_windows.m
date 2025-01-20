function test_windows()

    % Go to the directory of this repository.
    cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..'));

    solvers = {@fminsearch, @fminunc};
    options.solver_names = {'simplex', 'bfgs'};

    benchmark(solvers, options);
end