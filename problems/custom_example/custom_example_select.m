function problem_names = custom_example_select(options)
    % This is a toy example to show how to write a custom problem selector.
    % We will use the mat file `custom_example_info.mat` created by
    % `custom_example_getInfo` to select the problems.

    % Initialization.
    problem_names = {};

    % Set default values for options.
    if ~isfield(options, 'ptype')
        options.ptype = 'ubln';
    end
    if ~isfield(options, 'mindim')
        options.mindim = 1;
    end
    if ~isfield(options, 'maxdim')
        options.maxdim = Inf;
    end
    if ~isfield(options, 'minb')
        options.minb = 0;
    end
    if ~isfield(options, 'maxb')
        options.maxb = Inf;
    end
    if ~isfield(options, 'mincon')
        options.mincon = 0;
    end
    if ~isfield(options, 'maxcon')
        options.maxcon = Inf;
    end
    if ~isfield(options, 'excludelist')
        options.excludelist = {};
    end

    % Load the data from a .mat file.
    load('custom_example_info.mat', 'problem_info');

    % Loop through each problem and check the criteria.
    for i_problem = 2:size(problem_info, 1)
        problem_name = problem_info{i_problem, 1};
        ptype = problem_info{i_problem, 2};
        dim = problem_info{i_problem, 3};
        mb = problem_info{i_problem, 4};
        m_con = problem_info{i_problem, 5};

        % Check the criteria.
        if ismember(ptype, options.ptype) && ...
           (dim >= options.mindim && dim <= options.maxdim) && ...
           (mb >= options.minb && mb <= options.maxb) && ...
           (m_con >= options.mincon && m_con <= options.maxcon) && ...
           ~ismember(problem_name, options.excludelist)
            problem_names{end + 1} = problem_name;
        end
    end
end