function [x, fval] = uobyqa_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.solver = 'uobyqa';

    [x, fval] = pdfo(problem);

end