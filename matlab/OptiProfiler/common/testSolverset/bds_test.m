function [x, fval] = bds_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval)

    options.maxfun = max_eval;
    [x, fval] = bds(fun, x0, options);

end