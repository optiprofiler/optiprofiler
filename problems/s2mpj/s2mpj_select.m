function [problem_names, argins] = s2mpj_select(options)
%S2MPJ_SELECT selects the problems in S2MPJ that satisfy given criteria.
%
%   Users only need to use the following signature to call this function:
%
%   PROBLEM_NAMES = S2MPJ_SELECT(OPTIONS) returns the names of selected
%   problems from S2MPJ that satisfy the criteria in OPTIONS as a cell array
%   PROBLEM_NAMES. More details about S2MPJ can be found in the official
%   website: <https://github.com/GrattonToint/S2MPJ>.
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
%       - minlcon: the minimum number of linear constraints of the problems to
%         be selected. Default is 0.
%       - maxlcon: the maximum number of linear constraints of the problems to
%         be selected. Default is Inf.
%       - minnlcon: the minimum number of nonlinear constraints of the problems
%         to be selected. Default is 0.
%       - maxnlcon: the maximum number of nonlinear constraints of the problems
%         to be selected. Default is Inf.
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
%   Three things to note:
%
%       1. All the information about the problems can be found in a csv file
%          named 'probinfo_matlab.csv' in the same directory as this function.
%       2. The problem name may appear in the form of 'problem_name_dim_mcon'
%          where 'problem_name' is the name of the problem, 'dim' is the
%          dimension of the problem, and 'mcon' is the number of linear and
%          nonlinear constraints of the problem. This case only happens when
%          this problem can accept extra arguments to change the dimension or
%          the number of constraints. This information is stored in the
%          'probinfo_matlab.csv' file as the last few columns.
%       3. There is a file `config.txt` in the same directory as this function.
%          This file can be used to set the options `variable_size` and
%          `test_feasibility_problems`. Details about these two options can be
%          found in the comments in the `config.txt` file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set the `variable_size` %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set default options in the config file if not provided.
    variable_size = 'default';
    test_feasibility_problems = 0;
    % Check whether there is a file named 'config.txt' in the current directory.
    config_path = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
    if exist(config_path, 'file') == 2
        % Read the content of the file.
        try
            fid = fopen(config_path, 'r');
            % Find the line starting with 'variable_size=' and 'test_feasibility_problems='
            while ~feof(fid)
                line = fgetl(fid);
                if startsWith(line, 'variable_size=')
                    variable_size = strtrim(extractAfter(line, 'variable_size='));
                elseif startsWith(line, 'test_feasibility_problems=')
                    test_feasibility_problems = str2double(strtrim(extractAfter(line, 'test_feasibility_problems=')));
                end
            end
            fclose(fid);
        catch
            warning('Failed to read the file `config.txt`. Using default options.');
        end
    end
    % Check the validity of `variable_size` and `test_feasibility_problems`.
    if ~((ischar(variable_size) || isstring(variable_size)) && ismember(char(variable_size), {'default', 'min', 'max', 'all'}))
        error("Invalid `variable_size` in the file `config.txt`. Please set it to 'default', 'min', 'max', or 'all' (without quotes).");
    end
    if ~ismember(test_feasibility_problems, [0, 1, 2])
        error("Invalid `test_feasibility_problems` in the file `config.txt`. Please set it to 0, 1, or 2.");
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialization.
    problem_names = {};
    argins = {};

    % Check whether the options are valid.
    valid_fields = {'ptype', 'mindim', 'maxdim', 'minb', 'maxb', 'minlcon', 'maxlcon', 'minnlcon', 'maxnlcon', 'mincon', 'maxcon', 'oracle', 'excludelist'};
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
    if ~isfield(options, 'minlcon')
        options.minlcon = 0;
    end
    if ~isfield(options, 'maxlcon')
        options.maxlcon = Inf;
    end
    if ~isfield(options, 'minnlcon')
        options.minnlcon = 0;
    end
    if ~isfield(options, 'maxnlcon')
        options.maxnlcon = Inf;
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

    % Add some problems to excludelist. These problems have bugs or are not suitable for benchmarking.
    % 'DANWOODLS' and 'MISRA1CLS' have bugs at this moment (bugs in SIF files where bound constraints are not defined).
    % 'ROSSIMP1', 'ROSSIMP2', and 'ROSSIMP3' are duplicated with 'ROSENBR'.
    options.excludelist = unique([options.excludelist, {'DANWOODLS', 'MISRA1CLS', 'ROSSIMP1', 'ROSSIMP2', 'ROSSIMP3'}]);

    % Load the data from a .mat file.
    load('probinfo_matlab.mat', 'probinfo');

    for i_problem = 2:size(probinfo, 1)
        problem_name = probinfo{i_problem, 1};
        ptype = probinfo{i_problem, 2};
        dim = probinfo{i_problem, 4};
        mb = probinfo{i_problem, 5};
        mlcon = probinfo{i_problem, 9};
        mnlcon = probinfo{i_problem, 10};
        mcon = probinfo{i_problem, 8};
        is_feasibility = probinfo{i_problem, 18};
        argin = probinfo{i_problem, 25};
        dims = probinfo{i_problem, 26};
        mbs = probinfo{i_problem, 27};
        mlcons = probinfo{i_problem, 31};
        mnlcons = probinfo{i_problem, 32};
        mcons = probinfo{i_problem, 30};

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

        % Check if the problem is a feasibility problem.
        if test_feasibility_problems == 0
            if is_feasibility
                continue;
            end
        elseif test_feasibility_problems == 1
            if ~is_feasibility
                continue;
            end
        % elseif test_feasibility_problems == 2
        %     Do nothing, include all problems.
        end

        % If the default dimension and number of constraints satisfy the criteria, we add the problem.
        default_satisfy = dim >= options.mindim && dim <= options.maxdim && mb >= options.minb && mb <= options.maxb && mlcon >= options.minlcon && mlcon <= options.maxlcon && mnlcon >= options.minnlcon && mnlcon <= options.maxnlcon && mcon >= options.mincon && mcon <= options.maxcon;
        if default_satisfy && (isempty(dims) || (strcmp(variable_size, 'default')))
            problem_names{end + 1} = problem_name;
            argins{end + 1} = {};
        end
        if strcmp(variable_size, 'default')
            continue;
        end

        % If `variable_size` is  'min', 'max', or 'all', we need to check the variable sizes.
        if ~isempty(dims)
            mask = (dims >= options.mindim & dims <= options.maxdim) & (mbs >= options.minb & mbs <= options.maxb) & (mlcons >= options.minlcon & mlcons <= options.maxlcon) & (mnlcons >= options.minnlcon & mnlcons <= options.maxnlcon) & (mcons >= options.mincon & mcons <= options.maxcon) & (dims ~= dim | mcons ~= mcon);
            idxs = find(mask);
            if isempty(idxs)
                if default_satisfy
                    problem_names{end + 1} = problem_name;
                    argins{end + 1} = {};
                end
                continue;
            end
            if strcmp(variable_size, 'min')
                % Find the minimal dimension among all configurations that satisfy the options
                min_dim = min(dims(idxs));
                idxs_dim = idxs(dims(idxs) == min_dim);
                % Among those, find the one with the minimal number of constraints.
                min_mcon = min(mcons(idxs_dim));
                idxs = idxs_dim(mcons(idxs_dim) == min_mcon);
                idxs = idxs(1);
                % Compare with the default configuration.
                if default_satisfy && (dim < min_dim || (dim == min_dim && mcon < min_mcon))
                    problem_names{end + 1} = problem_name;
                    argins{end + 1} = {};
                    continue;
                end
            elseif strcmp(variable_size, 'max')
                % Find the maximal dimension among all configurations that satisfy the options
                max_dim = max(dims(idxs));
                idxs_dim = idxs(dims(idxs) == max_dim);
                % Among those, find the one with the maximal number of constraints.
                max_mcon = max(mcons(idxs_dim));
                idxs = idxs_dim(mcons(idxs_dim) == max_mcon);
                idxs = idxs(1);
                % Compare with the default configuration.
                if default_satisfy && (dim > max_dim || (dim == max_dim && mcon > max_mcon))
                    problem_names{end + 1} = problem_name;
                    argins{end + 1} = {};
                    continue;
                end
            elseif strcmp(variable_size, 'all')
                % Include all configurations that satisfy the options.
                % Do nothing, we will process all idxs.
            end
            for k = 1:numel(idxs)
                idx = idxs(k);
                if mcons(idx) == 0
                    problem_names{end + 1} = [problem_name, '_', num2str(dims(idx))];
                else
                    problem_names{end + 1} = [problem_name, '_', num2str(dims(idx)), '_', num2str(mcons(idx))];
                end
                if iscell(argin)
                    argins{end + 1} = {argin{1:end-1}, argin{end}(idx)};
                else
                    argins{end + 1} = argin(idx);
                end
            end
        end
    end
end