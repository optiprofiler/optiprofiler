function custom_example_getInfo()
    % This is a toy example to show one way of collecting information about
    % the problems in the custom_example problem library.
    % We will store the information in a MAT file so that `custom_example_select`
    % can use it to select problems.

    % Collect all the m file names in the directory `matlab_problems`.
    mydir = fileparts(mfilename('fullpath'));
    mfiles = dir(fullfile(mydir, 'matlab_problems', '*.m'));
    mfile_names = {mfiles.name};
    mfile_names = regexprep(mfile_names, '\.m$', ''); % Remove the .m extension

    % Initial the cell array to store the problem information.
    n_problems = length(mfile_names);
    problem_info = cell(n_problems + 1, 5);
    problem_info(1, :) = {'name', 'ptype', 'dim', 'mb', 'mcon'};

    % Loop through each m file and collect the information.
    for i_problem = 1:n_problems
        % Get the problem name.
        problem_name = mfile_names{i_problem};
        % Load the problem.
        p = custom_example_load(problem_name);
        % Record the problem information.
        problem_info{i_problem + 1, 1} = p.name;
        problem_info{i_problem + 1, 2} = p.ptype;
        problem_info{i_problem + 1, 3} = p.n;
        problem_info{i_problem + 1, 4} = sum(~isinf(-p.xl)) + sum(~isinf(p.xu));
        problem_info{i_problem + 1, 5} = p.m_linear_eq + p.m_linear_ub + p.m_nonlinear_eq + p.m_nonlinear_ub;
    end

    % Save the problem information to a MAT file.
    save(fullfile(mydir, 'custom_example_info.mat'), 'problem_info');
end