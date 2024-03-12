function [x, fval] = nloptnewuoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)

    opt = struct();
    opt.algorithm = NLOPT_LN_NEWUOA;
    opt.min_objective = fun;

    [x, fval, ~] = nlopt_optimize(opt, x0);

end