function problem = matcutest_load(problem_name)
%MATCUTEST_LOAD coverts the MatCUTEst problem name to a Problem class instance.
%
%   PROBLEM = MATCUTEST_LOAD(PROBLEM_NAME) returns a Problem class instance
%   PROBLEM that corresponds to the problem named PROBLEM_NAME in MatCUTEst.
%   More details about MatCUTEst can be found in the official website:
%   <https://github.com/matcutest>.
%
%   You may use the function `matcutest_select` to get the problem names you
%   want.
%
%   Note that MatCUTEst is only available in Linux.
%

    % Convert 'problem_name' to a char.
    problem_name = char(problem_name);

    % Use functions from MatCUTEst to get the problem structure.
    try
        pb = macup(problem_name);
    catch
        error("MATLAB:matcutest_load:errormacup", "Error occurred while using macup function on problem %s. Please check if MatCUTEst is installed correctly.", problem_name);
    end

    % Get the objective function.
    fun = @(x) getfun(pb, x);

    % Get the gradient of the objective function.
    grad = @(x) getgrad(pb, x);

    % Get the Hessian of the objective function.
    hess = @(x) gethess(pb, x);

    % Get the initial guess, lower bound, and upper bound.
    x0 = pb.x0;
    xl = pb.lb;
    xu = pb.ub;

    % Get the linear constraints.
    aeq = pb.Aeq;
    beq = pb.beq;
    aub = pb.Aineq;
    bub = pb.bineq;

    % Get the nonlinear constraints and their Jacobians.
    ceq = @(x) getceq(pb, x);
    cub = @(x) getcub(pb, x);
    jceq = @(x) getjceq(pb, x);
    jcub = @(x) getjcub(pb, x);

    problem = Problem(struct('name', problem_name, 'fun', fun, 'grad', grad, 'hess', hess, 'x0', x0, 'xl', xl, 'xu', xu, 'aeq', aeq, 'beq', beq, 'aub', aub, 'bub', bub, 'ceq', ceq, 'cub', cub, 'jceq', jceq, 'jcub', jcub));
    
end

function fx = getfun(pb, x)
    % Get the objective function value.
    try
        evalc('fx = pb.objective(x)');
    catch
        fx = NaN(0, 1);
    end
end

function gx = getgrad(pb, x)
    % Get the gradient of the objective function.
    try
        evalc('[~, gx] = pb.objective(x)');
    catch
        gx = NaN(0, 1);
    end
end

function hx = gethess(pb, x)
    % Get the Hessian of the objective function.
    try
        evalc('[~, ~, hx] = pb.objective(x)');
    catch
        hx = NaN(0, 1);
    end
end

function cubx = getcub(pb, x)
    % Get the nonlinear inequality constraints.
    try
        evalc('cubx = pb.nonlcon(x)');
    catch
        cubx = NaN(0, 1);
    end
end

function ceqx = getceq(pb, x)
    % Get the nonlinear equality constraints.
    try
        evalc('[~, ceqx] = pb.nonlcon(x)');
    catch
        ceqx = NaN(0, 1);
    end
end

function jcubx = getjcub(pb, x)
    % Get the Jacobian of the nonlinear inequality constraints.
    try
        evalc('[~, ~, jcubx] = pb.nonlcon(x)');
        jcubx = jcubx';
    catch
        jcubx = NaN(0, 1);
    end
end

function jceqx = getjceq(pb, x)
    % Get the Jacobian of the nonlinear equality constraints.
    try
        evalc('[~, ~, ~, jceqx] = pb.nonlcon(x)');
        jceqx = jceqx';
    catch
        jceqx = NaN(0, 1);
    end
end