function problem_names = selector(options)
    % Set default values for options
    if ~isfield(options, 'problem_type')
        options.problem_type = 'ubln';
    end
    if ~isfield(options, 'mindim')
        options.mindim = 1;
    end
    if ~isfield(options, 'maxdim')
        options.maxdim = Inf;
    end
    if ~isfield(options, 'mincon')
        options.mincon = 0;
    end
    if ~isfield(options, 'maxcon')
        options.maxcon = Inf;
    end

    % Load the data from a .mat file
    load('probinfo.mat', 'problem_data');

    % Initialize the selection mask as all false
    selection_mask = false(1, size(problem_data, 1) - 1);

    for i_problem = 2:size(problem_data, 1)
        problem_type = problem_data{i_problem, 2};
        dim = problem_data{i_problem, 3};
        m_con = problem_data{i_problem, 7};

        % Check problem_type and dimension criteria
        selection_mask(i_problem - 1) = (ismember(problem_type, options.problem_type) && ...
            dim >= options.mindim && dim <= options.maxdim && ...
            m_con >= options.mincon && m_con <= options.maxcon);
    end

    % Extract names based on the selection mask
    problem_names = problem_data(2:end, 1)';
    problem_names = problem_names(selection_mask);
    
end