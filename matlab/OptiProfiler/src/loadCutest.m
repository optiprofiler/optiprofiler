function problem = loadCutest(problem_name, problem_options)
    %{
    Load a CUTEst problem.

    Parameters
    ----------
    problem_name : str
        Name of the CUTEst problem to load.

    Returns
    -------
    Problem
        Loaded problem.

    Other Parameters
    ----------------
    mindim : int, optional
        Minimum number of variables.
    maxdim : int, optional
        Maximum number of variables.
    mincon : int, optional
        Minimum number of linear and nonlinear constraints.
    maxcon : int, optional
        Maximum number of linear and nonlinear constraints.

    Raises
    ------
    ProblemError
        If the problem cannot be loaded.
    %}

    cutest_problem = [];

    if nargin < 1
        error("Not enough input arguments.");
    elseif ~ischar(problem_name)
        error("The problem name must be a string.");
    else
        % TODO: this part needs to be rewritten
        try
            cutest_problem = macup(problem_name);
        catch ME
            warning(ME.identifier, "Failed to load CUTEst problem " + problem_name + ".");
        end

        if isempty(cutest_problem)
            error("Failed to load CUTEst problem " + problem_name + ".");
        end

        if nargin == 2 && ~isempty(problem_options)
            optionKeys = fieldnames(problem_options);
            validKeys = {enumeration('ProblemOptionKey').value};
            for i = 1:numel(optionKeys)
                key = optionKeys{i};
                value = problem_options.(key);
                if ~ismember(key, validKeys)
                    error("Unknown option: " + key + ".");
                end
                if isfloat(value) && mod(value, 1) == 0
                    problem_options.(key) = int32(value);
                end
                if ~isinteger(problem_options.(key)) || problem_options.(key) < 0
                    error("The argument " + key + " must be a nonnegative integer.");
                end
            end

            if isfield(problem_options, ProblemOptionKey.N_MIN.value) && isfield(problem_options, ProblemOptionKey.N_MAX.value) && problem_options.(ProblemOptionKey.N_MIN.value) > problem_options.(ProblemOptionKey.N_MAX.value)
                error("The argument " + ProblemOptionKey.N_MIN.value + " must be less than or equal to " + ProblemOptionKey.N_MAX.value + ".");
            end

            if isfield(problem_options, ProblemOptionKey.M_MIN.value) && isfield(problem_options, ProblemOptionKey.M_MAX.value) && problem_options.(ProblemOptionKey.M_MIN.value) > problem_options.(ProblemOptionKey.M_MAX.value)
                error("The argument " + ProblemOptionKey.M_MIN.value + " must be less than or equal to " + ProblemOptionKey.M_MAX.value + ".");
            end

            if ~isempty(cutest_problem) && ~isValid(cutest_problem, problem_options)
                error("CUTEst problem " + problem_name + " is invalid.");
            end
        end
    end

    problem = 1;

end



function is_valid = isValid(cutest_problem, problem_options)

    % TODO: why need vtype and how?

    is_valid = true;
    cutest_problem_dim = size(cutest_problem.x0, 1);

    % Check that the dimensions are within the specified range.
    if isfield(problem_options, ProblemOptionKey.N_MIN.value)
        is_valid = is_valid && cutest_problem_dim >= problem_options.(ProblemOptionKey.N_MIN.value);
    end
    if isfield(problem_options, ProblemOptionKey.N_MAX.value)
        is_valid = is_valid && cutest_problem_dim <= problem_options.(ProblemOptionKey.N_MAX.value);
    end
    if isfield(problem_options, ProblemOptionKey.M_MIN.value)
        is_valid = is_valid && cutest_problem_dim >= problem_options.(ProblemOptionKey.M_MIN.value);
    end
    if isfield(problem_options, ProblemOptionKey.M_MAX.value)
        is_valid = is_valid && cutest_problem_dim <= problem_options.(ProblemOptionKey.M_MAX.value);
    end
end