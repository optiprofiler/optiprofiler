function testOptiProfiler()
%TESTOPTIPROFILER Do some tests to check if OptiProfiler is installed correctly.

    warning('off', 'all');

    fprintf("\nRunning tests for OptiProfiler...\n\n");

    % Test if problem from S2MPJ is callable.
    try
        p = s2mpj_load('BEALE');
        assert(isequal(p.x0, [1; 1]));
        assert(isequal(p.name, 'BEALE'));
        assert(isnumeric(p.fun(p.x0)) && isreal(p.fun(p.x0)) && isscalar(p.fun(p.x0)));
        fprintf("OptiProfiler successfully loaded problem 'BEALE'.\n\n");
    catch
        error("OptiProfiler FAILED a test: Unable to load problem 'BEALE'.");
    end

    % Test if function `benchmark` is callable.
    path = pwd;
    path_out = fullfile(path, 'out');
    try
        options.problem = s2mpj_load('BEALE');
        evalc("benchmark({@fminsearch_test1, @fminsearch_test2}, options)");
        try
            rmdir(path_out, 's');
        catch
        end
        fprintf("OptiProfiler successfully ran function 'benchmark'.\n\n");
    catch
        try
            rmdir(path_out, 's');
        catch
        end
        error("OptiProfiler FAILED a test: Unable to call function 'benchmark'.");
    end

    % Failed checks throw above, so reaching here means both checks passed.
    fprintf("OptiProfiler passed all tests.\n\n");
    fprintf("You may now try 'help optiprofiler' for information on the usage of the package.\n\n");

    warning('on', 'all');
end

function x = fminsearch_test1(fun, x0)

    options = optimset('MaxFunEvals', 200);
    x = fminsearch(fun, x0, options);
end

function x = fminsearch_test2(fun, x0)

    options = optimset('MaxFunEvals', 500);
    x = fminsearch(fun, x0, options);
end
