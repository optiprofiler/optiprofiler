function [problem_names, argins] = s_select(options)
%S_SELECT selects the problems in S2MPJ that satisfy given criteria.
%
%   Users only need to use the following signature to call this function:
%
%   PROBLEM_NAMES = S_SELECT(OPTIONS) returns the names of selected problems
%   from S2MPJ that satisfy the criteria in OPTIONS as a cell array
%   PROBLEM_NAMES.
%
%   OPTIONS is a struct with the following fields:
%
%       - ptype: the type of the problems to be selected. It should be a string
%         or char consisting of any combination of 'u' (unconstrained), 'b'
%         (bound constrained), 'l' (linearly constrained), and 'n' (nonlinearly
%         constrained), such as 'b', 'ul', 'ubn'. Default is 'ubln'.
%       - mindim: the minimum dimension of the problems to be selected. Default
%         is 1.
%       - maxdim: the maximum dimension of the problems to be selected. Default
%         is Inf.
%       - minb: the minimum number of bound constraints of the problems to be
%         selected. Default is 0.
%       - maxb: the maximum number of bound constraints of the problems to be
%         selected. Default is Inf.
%       - mincon: the minimum number of linear and nonlinear constraints of the
%         problems to be selected. Default is 0.
%       - maxcon: the maximum number of linear and nonlinear constraints of the
%         problems to be selected. Default is Inf.
%       - oracle: the oracle provided by the problem. If it is 0, then the
%         problem should provide zeroth-order information. If it is 1, then the
%         problem should provide both zeroth-order and first-order information.
%         If it is 2, then the problem should provide zeroth-order,
%         first-order, and second-order information. Default is 0.
%       - excludelist: the list of problems to be excluded. Default is not to
%         exclude any problem.
%
%   Two things to note:
%       1. All the information about the problems can be found in a csv file
%          named 'probinfo.csv' in the same directory as this function.
%       2. The problem name may appear in the form of 'problem_name_dim_m_con'
%          where 'problem_name' is the name of the problem, 'dim' is the
%          dimension of the problem, and 'm_con' is the number of linear and
%          nonlinear constraints of the problem. This case only happens when
%          this problem can accept extra arguments to change the dimension or
%          the number of constraints. This information is stored in the
%          'probinfo.csv' file as the last few columns.
%

    % Initialization.
    problem_names = {};
    argins = {};

    % Check whether the options are valid.
    valid_fields = {'ptype', 'mindim', 'maxdim', 'minb', 'maxb', 'mincon', 'maxcon', 'oracle', 'excludelist'};
    if ~isstruct(options) || (~isempty(fieldnames(options)) && ~all(ismember(fieldnames(options), valid_fields)))
        error('The input argument `options` is invalid.');
    end

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
    if ~isfield(options, 'oracle')
        options.oracle = 0;
    end
    if ~isfield(options, 'excludelist')
        options.excludelist = {};
    end

    % Load the data from a .mat file.
    load('probinfo.mat', 'probinfo');

    for i_problem = 2:size(probinfo, 1)
        problem_name = probinfo{i_problem, 1};
        ptype = probinfo{i_problem, 2};
        dim = probinfo{i_problem, 4};
        mb = probinfo{i_problem, 5};
        m_con = probinfo{i_problem, 8};
        argin = probinfo{i_problem, 25};
        dims = probinfo{i_problem, 26};
        mbs = probinfo{i_problem, 27};
        m_cons = probinfo{i_problem, 30};

        % If the oracle is not 0, then we exclude problem 'NOZZLEfp' since it does not have first- or second-order information.
        % "NOZZLEfp" is designed to simulate jet impingement cooling.
        % See https://optimization-online.org/wp-content/uploads/2024/03/Design_Optimization_Of_A_Jet_Plate_for_Impingement_Cooling-1.pdf
        if options.oracle ~= 0
            options.excludelist = [options.excludelist, 'NOZZLEfp'];
        end

        % Check if the problem is in the exclude list.
        if ~isempty(options.excludelist) && ismember(problem_name, options.excludelist)
            continue;
        end

        % Check if the problem type satisfies the criteria.
        if ~ismember(ptype, options.ptype)
            continue;
        end

        % If the default dimension and number of constraints satisfy the criteria, we add the problem.
        % That means, we will not consider the changeable dimension and number of constraints if the
        % default dimension and number of constraints satisfy the criteria.
        if dim >= options.mindim && dim <= options.maxdim && mb >= options.minb && mb <= options.maxb && m_con >= options.mincon && m_con <= options.maxcon
            problem_names{end + 1} = problem_name;
            argins{end + 1} = {};
            continue;
        end

        % If the default dimension and number of constraints do not satisfy the criteria, we consider
        % the changeable dimension and number of constraints.
        if ~isempty(dims) && any(dims >= options.mindim & dims <= options.maxdim) && any(mbs >= options.minb & mbs <= options.maxb) && any(m_cons >= options.mincon & m_cons <= options.maxcon)
            idx = find(dims >= options.mindim & dims <= options.maxdim, 1, 'first');
            if m_cons(idx) == 0
                problem_names{end + 1} = [problem_name, '_', num2str(dims(idx))];
            else
                problem_names{end + 1} = [problem_name, '_', num2str(dims(idx)), '_', num2str(m_cons(idx))];
            end
            if iscell(argin)
                argins{end + 1} = {argin{1:end-1}, argin{end}(idx)};
            else
                argins{end + 1} = argin(idx);
            end
        end
    end
end