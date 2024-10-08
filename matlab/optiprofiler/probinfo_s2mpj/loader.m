function problem = loader(problem_name)
    
    funcHandle = str2func(problem_name);
    fun = @(x) feval(@(y) y{1}, num2cell(funcHandle('fx', x)));

    p = funcHandle('setup');
    x0 = p.x0;
    xl = p.xlower;
    xu = p.xupper;
    
    if ~isfield(p, 'lincons')
        p.lincons = [];
    end
    if ~isfield(p, 'nle')
        p.nle = 0;
    end
    if ~isfield(p, 'neq')
        p.neq = 0;
    end
    if ~isfield(p, 'nge')
        p.nge = 0;
    end

    p.nonlincons = setdiff(1:p.m, p.lincons);

    try
        [~] = evalc('[cx, Jx] = funcHandle("cJx", x0)');
        bx = Jx * x0 - cx;
        constraint = @(x) funcHandle('cx', x);
    catch
        Jx = NaN(0, size(x0, 1));
        bx = NaN(0, 1);
        constraint = @(x) NaN(0, 1);
    end

    idx_le = 1:p.nle;
    idx_eq = p.nle + 1 : p.nle + p.neq;
    idx_ge = p.nle + p.neq + 1 : p.nle + p.neq + p.nge;
    idx_aeq = intersect(idx_eq, p.lincons);
    idx_aub_le = intersect(idx_le, p.lincons);
    idx_aub_ge = intersect(idx_ge, p.lincons);
    idx_cle = intersect(idx_le, p.nonlincons);
    idx_ceq = intersect(idx_eq, p.nonlincons);
    idx_cge = intersect(idx_ge, p.nonlincons);
    aeq = Jx(idx_aeq, :);
    aub = [Jx(idx_aub_le, :); -Jx(idx_aub_ge, :)];
    beq = bx(idx_aeq);
    bub = [bx(idx_aub_le); -bx(idx_aub_ge)];
    m_nonlinear_ub = length(idx_cle) + length(idx_cge);
    m_nonlinear_eq = length(idx_ceq);

    getidx = @(x, idx) x(idx);
    cub = @(x) [getidx(constraint(x), idx_cle), -getidx(constraint(x), idx_cge)];
    ceq = @(x) getidx(constraint(x), idx_ceq);

    % Debug S2MPJ: unconstrained problems are added bound constraints (x>=0) by mistake
    
    % Check if p.m == 0 & xl = zeros(size(x0)) and xu = Inf(size(x0))
    if p.m == 0 && all(xl == 0) && all(xu == Inf)
        xl = -Inf(size(x0));
    end
    
    problem = Problem(struct('name', problem_name, 'fun', fun, 'x0', x0, 'xl', xl, 'xu', xu, 'aeq', aeq, 'beq', beq, 'aub', aub, 'bub', bub, 'cub', cub, 'ceq', ceq, 'm_nonlinear_ub', m_nonlinear_ub, 'm_nonlinear_eq', m_nonlinear_eq));
    
end