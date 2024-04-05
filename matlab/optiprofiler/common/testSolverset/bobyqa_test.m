function [x, fval] = bobyqa_test(fun, x0, xl, xu)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.lb = xl;
    problem.ub = xu;
    problem.solver = 'bobyqa';

    [x, fval] = pdfo(problem);
end