function [x, fval] = bobyqa_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.lb = xl;
    problem.ub = xu;
    problem.options.maxfun = max_eval;
    problem.solver = 'bobyqa';

    [x, fval] = pdfo(problem);
end