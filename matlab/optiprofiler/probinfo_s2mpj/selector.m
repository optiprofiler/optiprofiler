function problem_names = selector(options)
    % Set default values for options
    if ~isfield(options, 'type')
        options.type = 'u';
    end
    if ~isfield(options, 'mindim')
        options.mindim = 1;
    end
    if ~isfield(options, 'maxdim')
        options.maxdim = 5;
    end
    if ~isfield(options, 'mincon')
        options.mincon = 0;
    end
    if ~isfield(options, 'maxcon')
        options.maxcon = 0;
    end

    % Load the data from a .mat file
    load('probinfo.mat', 'problem_data');

    % Initialize the selection mask as all false
    selection_mask = false(1, size(problem_data, 1) - 1);

    for i_problem = 2:size(problem_data, 1)
        type = problem_data{i_problem, 2};
        dim = problem_data{i_problem, 3};
        m = problem_data{i_problem, 4};

        % Check type and dimension criteria
        if ismember(type, options.type) && dim >= options.mindim && dim <= options.maxdim && m >= options.mincon && m <= options.maxcon
            selection_mask(i_problem - 1) = true;
        end
    end

    % Extract names based on the selection mask
    problem_names = problem_data(2:end, 1)';
    problem_names = problem_names(selection_mask);
    
end