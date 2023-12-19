function [x, fval] = fminsearch_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval)

    n = length(x0);
    options = optimset('MaxFunEvals', max_eval);
    [x, fval] = fminsearch(fun, x0, options);

end