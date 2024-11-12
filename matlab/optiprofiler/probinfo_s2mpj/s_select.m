function [problem_names, argins] = s_select(options)
%S_SELECT specific problem selector for the problem set "S2MPJ".

    % Initialization
    problem_names = {};
    argins = {};

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
    load('probinfo.mat', 'probinfo');

    for i_problem = 2:size(probinfo, 1)
        problem_name = probinfo{i_problem, 1};
        problem_type = probinfo{i_problem, 2};
        dim = probinfo{i_problem, 4};
        m_con = probinfo{i_problem, 8};
        argin = probinfo{i_problem, 17};
        dims = probinfo{i_problem, 18};

        if ~ismember(problem_type, options.problem_type)
            continue;
        end
        if m_con < options.mincon || m_con > options.maxcon
            continue;
        end
        % If the default dimension satisfies the criteria, we add the problem
        if dim >= options.mindim && dim <= options.maxdim
            problem_names{end + 1} = problem_name;
            argins{end + 1} = {};
            continue;
        end
        % If the default dimension does not satisfy the criteria but there exists a changeable dimension that satisfies the criteria, we add the problem and the corresponding argin that leads to the smallest dimension that satisfies the criteria
        if ~isempty(dims) && any(dims >= options.mindim & dims <= options.maxdim)
            idx = find(dims >= options.mindim & dims <= options.maxdim, 1, 'first');
            problem_names{end + 1} = [problem_name, '_', num2str(dims(idx))];
            if iscell(argin)
                argins{end + 1} = {argin{1:end-1}, argin{end}(idx)};
            else
                argins{end + 1} = argin(idx);
            end
        end
    end
    
end