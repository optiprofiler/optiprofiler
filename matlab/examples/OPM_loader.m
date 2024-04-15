function problem = OPM_loader(problem_name)
    
    funcHandle = str2func(problem_name);
    fun = @(x) feval(@(y) y{1}, num2cell(funcHandle('objf', x)));
    [x0, ~, ~, xl, xu] = funcHandle('setup');

    if isempty(xl)
        if isempty(xu)
            problem = Problem(struct('fun', fun, 'x0', x0));
        else
            problem = Problem(struct('fun', fun, 'x0', x0, 'xu', xu));
        end
    else
        if isempty(xu)
            problem = Problem(struct('fun', fun, 'x0', x0, 'xl', xl));
        else
            problem = Problem(struct('fun', fun, 'x0', x0, 'xl', xl, 'xu', xu));
        end
    end
    
end