function [x, fval] = fminsearch_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval)

    n = length(x0);
    options = optimset('MaxFunEvals', floor(max_eval/n));
    [x, fval] = fminsearch(fun, x0, options);

end