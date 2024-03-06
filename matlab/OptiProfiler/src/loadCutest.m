function problem = loadCutest(problem_name, problem_options)
%LOADCUTEST loads a CUTEst problem.

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
            fprintf("Loading CUTEst problem " + problem_name + ".\n");
            cutest_problem = macup(problem_name);
        catch ME
            warning(ME.identifier, "Failed to load CUTEst problem " + problem_name + ".");
        end

        if isempty(cutest_problem)
            error("Failed to load CUTEst problem " + problem_name + ".");
        end

        problem_struct = struct();

        problem_struct.fun = @(x) cutest_problem.objective(x);
        problem_struct.x0 = cutest_problem.x0;
        
        xl = cutest_problem.lb;
        xl(xl <= -1e20) = -Inf;
        xu = cutest_problem.ub;
        xu(xu >= 1e20) = Inf;

        problem_struct.xl = xl;
        problem_struct.xu = xu;

        if cutest_problem.numcon > 0
            if cutest_problem.numlcon > 0
                if cutest_problem.numleq > 0
                    problem_struct.aeq = cutest_problem.Aeq;
                    problem_struct.beq = cutest_problem.beq;
                elseif cutest_problem.numlineq > 0
                    problem_struct.aub = cutest_problem.Aineq;
                    problem_struct.bub = cutest_problem.bineq;
                end
            elseif cutest_problem.numnlcon > 0
                if cutest_problem.numnleq > 0
                    problem_struct.ceq = @(x) selectNonlinear(cutest_problem, x, 1);
                    problem_struct.m_nonlinear_eq = cutest_problem.numnleq;
                elseif cutest_problem.numnlineq > 0
                    problem_struct.cub = @(x) selectNonlinear(cutest_problem, x, 0);
                    problem_struct.m_nonlinear_ub = cutest_problem.numnlineq;
                end
            end
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
                if ~isnumeric(problem_options.(key)) || ~(mod(value, 1) == 0) || problem_options.(key) < 0
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
                error("CUTEst problem " + problem_name + " is successfully loaded but invalid; it is discarded.");
            end

        end
    end

    problem = Problem(problem_struct);

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


function nonlinearVal = selectNonlinear(cutest_problem, x, number)

    % Since the function "nonlcon" in matcutest has four outputs ([nlcineq0, nlceq0, gnlcineq0, gnlceq0] = nonlcon(x0)), we need to specify nonlinear inequality and equality constraints separately.
    % number = 0: nonlinear inequality constraints
    % number = 1: nonlinear equality constraints

    if number == 0
        nonlinearVal = cutest_problem.nonlcon(x);
    elseif number == 1
        [~, nonlinearVal] = cutest_problem.nonlcon(x);
    else
        error("The argument number must be 0 or 1.");
    end
end
