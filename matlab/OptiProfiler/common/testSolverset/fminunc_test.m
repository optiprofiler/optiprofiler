function [x, fval] = fminunc_test(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval)

    n = length(x0);
    options = optimset('MaxFunEvals', floor(max_eval/n));
    [x, fval] = fminunc(fun, x0, options);

end