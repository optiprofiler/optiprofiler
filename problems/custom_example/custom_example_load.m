function problem = custom_example_load(problem_name)
    % This is a toy example to show how to write a custom problem loader.

    % Add the path to the problem directory.
    mydir = fileparts(mfilename('fullpath'));
    addpath(fullfile(mydir, 'matlab_problems'));

    % Load the problem.
    switch problem_name
        case 'custom1'
            problem = custom1();
        case 'custom2'
            problem = custom2();
        otherwise
            error('Unknown problem name: %s', problem_name);
    end
end