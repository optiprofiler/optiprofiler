function problem_names = matcutest_select(options)
%MATCUTEST_SELECT selects the problems in MatCUTEst that satisfy given criteria.
%
%   Users only need to use the following signature to call this function:
%
%   PROBLEM_NAMES = MATCUTEST_SELECT(OPTIONS) returns the names of selected
%   problems from MatCUTEst that satisfy the criteria in OPTIONS as a cell
%   array PROBLEM_NAMES. More details about MatCUTEst can be found in the
%   official website: <https://github.com/matcutest>.
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
%       - excludelist: the list of problems to be excluded. Default is not to
%         exclude any problem.
%
%   Two things to note:
%
%       1. MatCUTEst is only available in Linux.
%       2. There is a file `config.txt` in the same directory as this function.
%          This file can be used to set the option `test_feasibility_problems`.
%          Details about this option can be found in the comments in the
%          `config.txt` file.

    % Check whether the options are valid.
    valid_fields = {'ptype', 'mindim', 'maxdim', 'minb', 'maxb', 'minlcon', 'maxlcon', 'minnlcon', 'maxnlcon', 'mincon', 'maxcon', 'excludelist'};
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
    if ~isfield(options, 'excludelist')
        options.excludelist = {};
    end

    options.type = options.ptype;
    options = rmfield(options, 'ptype');
    options.blacklist = options.excludelist;
    options = rmfield(options, 'excludelist');

    % Read the configuration file to set options.
    test_feasibility_problems = 0;
    config_path = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
    if exist(config_path, 'file') == 2
        % Read the content of the file.
        try
            fid = fopen(config_path, 'r');
            % Find the line starting with 'test_feasibility_problems='
            while ~feof(fid)
                line = fgetl(fid);
                if startsWith(line, 'test_feasibility_problems=')
                    test_feasibility_problems = str2double(strtrim(extractAfter(line, 'test_feasibility_problems=')));
                end
            end
            fclose(fid);
        catch
            warning('Failed to read the file `config.txt`. Using default options.');
        end
    end
    % Check the validity of `test_feasibility_problems`.
    if ~ismember(test_feasibility_problems, [0, 1, 2])
        error("Invalid `test_feasibility_problems` in the file `config.txt`. Please set it to 0, 1, or 2.");
    end

    % Set the option `is_feasibility` depending on `test_feasibility_problems`.
    if test_feasibility_problems == 0
        options.is_feasibility = false;
    elseif test_feasibility_problems == 1
        options.is_feasibility = true;
    % Else do nothing.
    end

    % Use functions from MatCUTEst to select the problems.
    try
        problem_names = secup(options);
    catch
        error("MATLAB:matcutest_select:errorsecup", "Error occurred while using secup function. Please check if MatCUTEst is installed correctly.");
    end
end