function problem_names = custom_example_select(options)
    % This is a toy example to show how to write a custom problem selector.

    % Set default values for options.
    if ~isfield(options, 'ptype')
        options.ptype = 'u';
    end

    % Select the problems based on the options.
    problem_names = {};
    if ismember('u', options.ptype)
        problem_names = [problem_names, 'custom1'];
    end
    if ismember('b', options.ptype)
        problem_names = [problem_names, 'custom2'];
    end
end