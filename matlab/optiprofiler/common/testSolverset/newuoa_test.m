function [x, fval] = newuoa_test(fun, x0)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.solver = 'newuoa';

    [x, fval] = pdfo(problem);

end