function custom_getInfo()
    % This is a toy example to show one way of collecting information about
    % the problems in the custom problem library.
    % We will store the information in a MAT file so that `custom_select`
    % can use it to select problems.

    % Collect all the m file names in the directory `matlab_problems`.
    mydir = fileparts(mfilename('fullpath'));
    mfiles = dir(fullfile(mydir, 'matlab_problems', '*.m'));
    mfile_names = {mfiles.name};
    mfile_names = regexprep(mfile_names, '\.m$', ''); % Remove the .m extension
    mfile_names = sort(mfile_names); % Sort the file names.

    % Initial the cell array to store the problem information.
    n_problems = length(mfile_names);
    probinfo = cell(n_problems + 1, 7);
    probinfo(1, :) = {'name', 'ptype', 'dim', 'mb', 'mlcon', 'mnlcon', 'mcon'};

    % Loop through each m file and collect the information.
    for i_problem = 1:n_problems
        % Get the problem name.
        problem_name = mfile_names{i_problem};
        % Load the problem.
        p = custom_load(problem_name);
        % Record the problem information.
        probinfo{i_problem + 1, 1} = p.name;
        probinfo{i_problem + 1, 2} = p.ptype;
        probinfo{i_problem + 1, 3} = p.n;
        probinfo{i_problem + 1, 4} = p.mb;
        probinfo{i_problem + 1, 5} = p.mlcon;
        probinfo{i_problem + 1, 6} = p.mnlcon;
        probinfo{i_problem + 1, 7} = p.mcon;
    end

    % Save the problem information to a MAT file.
    save(fullfile(mydir, 'probinfo_matlab.mat'), 'probinfo');
end