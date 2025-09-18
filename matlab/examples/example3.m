function example3()
%EXAMPLE3 corresponds to "Example 3: useful option load" in the "Usage for MATLAB" part of our
% official website (www.optprof.com).
%
% This example shows how to use the `load` option in OptiProfiler.

    % Before running this example, you need to run `example2` successfully.
    % When the job is done, you can try this example.
    % The option `load` allows you to load the results from a previous benchmark run.

    % Print the information about this example.
    fprintf('\nThis is an example to load the results from a previous benchmark run (exmaple 2 if you just run it) and create the profiles on a selected problem set.\n');
    pause(1.5);
    fprintf('\nStart Example 3...\n\n');

    % Start example 3.
    options.load = 'latest';    % Load results from the latest run. You may set it to a specific timestamp (e.g., '20250528_163340') which is the last 15 characters of the folder name saving the results.
    options.mindim = 3;         % Select problems with dimension at least 3 from the loaded results.
    options.maxdim = 4;         % Select problems with dimension at most 4 from the loaded results.
    options.solvers_to_load = [1, 2];    % Load results of the first two solvers.

    scores = benchmark(options)
end