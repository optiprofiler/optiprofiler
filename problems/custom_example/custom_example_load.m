function problem = custom_example_load(problem_name)
    % This is a toy example to show how to write a custom problem loader.

    % Add the path to the problem directory.
    mydir = fileparts(mfilename('fullpath'));
    addpath(fullfile(mydir, 'matlab_problems'));

    % Load the problem by directly calling the problem's m file.
    function_handle = str2func(problem_name);
    problem = function_handle();
end