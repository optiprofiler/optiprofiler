function problem = s_load(problem_name, varargin)
%S_LOAD specific problem loader for the problem set "S2MPJ".
%   Problems in S2MPJ are in the following form.
%       min f(x)
%       s.t. xlower <= x <= xupper  (xlower may be -Inf and xupper may be Inf)
%            ci(x) <= 0, i = 1, ..., nle,
%            ci(x) = 0, i = nle + 1, ..., nle + neq,
%            ci(x) >= 0, i = nle + neq + 1, ..., nle + neq + nge.
%   We call all ci(x) as general constraints, which can be linear or nonlinear.
%
%   The problem is loaded by 'setup' action in the problem function.
%       pb = problem_name('setup');
%   The problem structure pb contains the following main fields which we will
%   use in this function.
%       pb.n: dimension of the problem.
%       pb.nle: number of less-than-or-equal-to constraints.
%       pb.neq: number of equality constraints.
%       pb.nge: number of greater-than-or-equal-to constraints.
%       pb.m: total number of constraints (m = nle + neq + nge).
%       pb.lincons: indices of linear constraints in {1, ..., m}.
%       pb.x0: initial guess.
%       pb.xlower: lower bound.
%       pb.xupper: upper bound.
%       pb.xtype: type of variables ('r' (real), 'i'(integer), or 'b'(binary)).
%
%   The objective function will be evaluated by 'fx' action.
%       fx = problem_name('fx', x);
%   The general constraints and Jacobians will be evaluated by 'cJx' action.
%       [cx, Jx] = problem_name('cJx', x);
%   Note that we need to use Jacobians Jx to get the linear constraints.
%
%   We will use all the above information and functions to create a instance of
%   the Problem class containing the following fields.
%       name: name of the problem, same as problem_name. If there are varargin,
%             then it means that the problem is parameterized. Thus, in this
%             case, the name will be 'problem_name_n_m', where n is the
%             dimension of the problem and m is the number of constraints. When
%             m =0, the name will be 'problem_name_n'.
%       fun: function handle to evaluate the objective function, which should
%            be a function handle calling the 'fx' action.
%       x_type: type of variables, same as pb.xtype.
%       x0: initial guess, same as pb.x0.
%       xl: lower bound, same as pb.xlower.
%       xu: upper bound, same as pb.xupper.
%       aeq: matrix of equality constraint coefficients, which should be
%            Jx(idx_aeq, :) with idx_aeq = intersect(idx_eq, lincons) and
%            idx_eq = nle + 1 : nle + neq.
%       beq: vector of equality constraint values, which should be bx(idx_aeq)
%            with bx = Jx * x - cx (since we want aeq @ x = beq).
%       aub: matrix of inequality constraint coefficients, which should be
%            [Jx(idx_aub_le, :); -Jx(idx_aub_ge, :)] with
%            idx_aub_le = intersect(idx_le, lincons), idx_le = 1 : nle,
%            idx_aub_ge = intersect(idx_ge, lincons), and
%            idx_ge = nle + neq + 1 : nle + neq + nge.
%       bub: vector of inequality constraint values, which should be
%            [bx(idx_aub_le); -bx(idx_aub_ge)] (since we want aub @ x <= bub).
%       ceq: function handle to evaluate the nonlinear equality constraints,
%            the value of which at x should be cx(idx_ceq) with
%            idx_ceq = intersect(idx_eq, nonlincons) and
%            nonlincons = {1, ..., m} \ lincons.
%       cub: function handle to evaluate the nonlinear inequality constraints,
%            the value of which at x should be [cx(idx_cle), -cx(idx_cge)] with
%            idx_cle = intersect(idx_le, nonlincons),
%            idx_cge = intersect(idx_ge, nonlincons).
%
%   The function will be used in two parts. One appears in 's_getInfo.m' to
%   load all the problems to collect the information. The other appears in
%   'solveOneProblem.m' (MATLAB source file of OptiProfiler) and plays the
%   role of default problem set loader.
%
%   In 's_getInfo.m', 's_load' may receive varargin to load the problem when
%   the problem is parameterized. However, in 'solveOneProblem.m', this will
%   not happen since the selector 's_select' will first convert the names of 
%   parameterized problems to "problem_name_n" and then call 's_load' without
%   varargin.

    % Convert 'problem_name' to a char.
    problem_name = char(problem_name);

    % Check if 'problem_name' has the pattern '_n_m' or '_n_0'.
    [is_problem_parameterized, problem_name, dim, m_con] = isproblem_parameterized(problem_name);

    if is_problem_parameterized
        load('probinfo.mat', 'probinfo');
        idx_pb = find(strcmp(probinfo(:, 1), problem_name), 1, 'first');
        if isempty(idx_pb)
            error('Problem %s not found in probinfo.mat.', problem_name);
        end
        % Find the corresponding argins for the specific dimension.
        dims = probinfo{idx_pb, 19};
        m_cons = probinfo{idx_pb, 20};
        argins = probinfo{idx_pb, 18};
        idx_dim = find(dims == dim & m_cons == m_con, 1, 'first');
        if iscell(argins)
            argin = {argins{1:end-1}, argins{end}(idx_dim)};
        else
            argin = argins(idx_dim);
        end

        % Load the problem with the specific dimension.
        if iscell(argin)
            problem = s_load(problem_name, argin{:});
        else
            problem = s_load(problem_name, argin);
        end
        return;
    end

    % % Check whether there exists 'problem_name.m' in some folder called 'matlab_problems'.
    % current_path = fileparts(mfilename('fullpath'));
    % problem_path = fullfile(current_path, 'matlab_problems', [problem_name, '.m']);
    % if ~exist(problem_path, 'file')
    %     error('Problem %s not found in the directory %s.', problem_name, problem_path);
    % end

    % % Specify which 'problem_name.m' to load.
    % addpath(fullfile(current_path, 'matlab_problems'));

    % Convert the problem name to a function handle for later use.
    funcHandle = str2func(problem_name);

    % Get the structure of the problem information.
    pb = funcHandle('setup', varargin{:});
    % problem_name = pb.name;

    if nargin < 2
        name = problem_name;
    else
        name = [problem_name, '_', num2str(pb.n)];
    end

    % Get the objective function.
    fun = @(x) getfun(problem_name, x);

    % Get the gradient of the objective function.
    grad = @(x) getgrad(problem_name, x);

    % Get the Hessian of the objective function.
    hess = @(x) gethess(problem_name, x);

    % Get the initial guess, lower bound, and upper bound.
    if isfield(pb, 'xtype')
        switch pb.xtype
            case 'r'
                x_type = 'real';
            case 'i'
                x_type = 'integer';
            case 'b'
                x_type = 'binary';
            otherwise
                x_type = 'real';
        end
    else
        x_type = 'real';
    end
    x0 = pb.x0;
    xl = pb.xlower;
    xu = pb.xupper;
    
    if ~isfield(pb, 'lincons')
        pb.lincons = [];
    end
    if ~isfield(pb, 'nle')
        pb.nle = 0;
    end
    if ~isfield(pb, 'neq')
        pb.neq = 0;
    end
    if ~isfield(pb, 'nge')
        pb.nge = 0;
    end

    try
        [~] = evalc('[cx, Jx] = funcHandle("cJx", x0)');
        bx = Jx * x0 - cx;
    catch
        Jx = NaN(0, size(x0, 1));
        bx = NaN(0, 1);
    end

    nonlincons = setdiff(1:pb.m, pb.lincons);
    idx_le = 1:pb.nle;
    idx_eq = pb.nle + 1 : pb.nle + pb.neq;
    idx_ge = pb.nle + pb.neq + 1 : pb.nle + pb.neq + pb.nge;
    idx_aeq = intersect(idx_eq, pb.lincons);
    idx_aub_le = intersect(idx_le, pb.lincons);
    idx_aub_ge = intersect(idx_ge, pb.lincons);
    idx_cle = intersect(idx_le, nonlincons);
    idx_ceq = intersect(idx_eq, nonlincons);
    idx_cge = intersect(idx_ge, nonlincons);
    aeq = Jx(idx_aeq, :);
    aub = [Jx(idx_aub_le, :); -Jx(idx_aub_ge, :)];
    beq = bx(idx_aeq);
    bub = [bx(idx_aub_le); -bx(idx_aub_ge)];

    getidx = @(y, idx) y(idx);
    ceq = @(x) getidx(getcx(problem_name, x), idx_ceq);
    cub = @(x) [getidx(getcx(problem_name, x), idx_cle); -getidx(getcx(problem_name, x), idx_cge)];
    Hceq = @(x) getidx(getHx(problem_name, x), idx_ceq);
    Hcub = @(x) [getidx(getHx(problem_name, x), idx_cle), getidx(getHx(problem_name, x), idx_cge)];
    
    getidx_mat = @(y, idx) y(idx, :);
    Jceq = @(x) getidx_mat(getJx(problem_name, x), idx_ceq);
    Jcub = @(x) [getidx_mat(getJx(problem_name, x), idx_cle); -getidx_mat(getJx(problem_name, x), idx_cge)];
    
    problem = Problem(struct('name', name, 'fun', fun, 'grad', grad, 'hess', hess, 'x_type', x_type, 'x0', x0, 'xl', xl, 'xu', xu, 'aeq', aeq, 'beq', beq, 'aub', aub, 'bub', bub, 'ceq', ceq, 'cub', cub, 'Jceq', Jceq, 'Jcub', Jcub, 'Hceq', Hceq, 'Hcub', Hcub));
    
end

function fx = getfun(problem_name, x)
    
    funcHandle = str2func(problem_name);
    try
        evalc("fx = funcHandle('fx', x)");
    catch
        % If the function evaluation fails, return NaN.
        % Note that some problems are feasibility problems, which do not have an objective function.
        % For these problems, we return 0.
        fx = 0;
    end
end

function gx = getgrad(problem_name, x)
    
    funcHandle = str2func(problem_name);
    try
        evalc("[~, gx] = funcHandle('fgx', x)");
        gx = full(gx);
    catch
        gx = NaN(0, 1);
    end
end

function Hx = gethess(problem_name, x)
    
    funcHandle = str2func(problem_name);
    try
        evalc("[~, ~, Hx] = funcHandle('fgHx', x)");
        Hx = full(Hx);
    catch
        Hx = NaN(0, 1);
    end
end

function cx = getcx(problem_name, x)
    
    funcHandle = str2func(problem_name);
    try
        evalc("cx = funcHandle('cx', x)");
        cx = full(cx);
        % Check if cx is a row vector, if yes, transpose it to a column vector.
        if size(cx, 1) == 1
            cx = cx';
        end
    catch
        cx = NaN(0, 1);
    end
end

function Jx = getJx(problem_name, x)
    
    funcHandle = str2func(problem_name);
    try
        evalc("[~, Jx] = funcHandle('cJx', x)");
        Jx = full(Jx);
    catch
        Jx = NaN(0, 1);
    end
end

function Hx = getHx(problem_name, x)
    
    funcHandle = str2func(problem_name);
    try
        evalc("[~, ~, Hx] = funcHandle('cJHx', x)");
        for i = 1:size(Hx, 2)
            Hx{i} = full(Hx{i});
        end
    catch
        Hx = NaN(0, 1);
    end
end

function [result, problem_name, dim, m_con] = isproblem_parameterized(problem_name)
    % Check if 'problem_name' has the pattern '_n_m' or 'n'. If it has, find
    % the position of the pattern and return the dimension 'n' and the number
    % of constraints 'm'.

    % Find the position of the pattern '_n_m' or '_n'.
    [match, idx] = regexp(problem_name, '_\d+_\d+', 'match', 'start');
    if isempty(match)
        [match, idx] = regexp(problem_name, '_\d+', 'match', 'start');
        if isempty(match)
            result = false;
            dim = 0;
            m_con = 0;
        else
            problem_name = problem_name(1:idx-1);
            results = sscanf(match{1}, '_%d');
            dim = results(1);
            m_con = 0;
            result = true;
        end
        return;
    end

    % Find the dimension 'n' and the number of constraints 'm'.
    problem_name = problem_name(1:idx-1);
    results = sscanf(match{1}, '_%d_%d');
    dim = results(1);
    m_con = results(2);
    result = true;
end