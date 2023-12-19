function [x, fval] = newuoa_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval)

    options.maxfun = max_eval;
    [x, fval] = pdfo(fun, x0, options);

end