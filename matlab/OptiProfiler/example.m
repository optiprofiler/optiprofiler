function example()

    addpath("src/")

    problem_names = secup(struct('type', 'u', 'maxdim', 5, 'mindim', 1));
    solvers = {@newuoa_test, @uobyqa_test};
    labels = {"newuoa", 'uobyqa'};

    runBenchmark(solvers, labels, problem_names, 'plain')

end

function [x, fval] = newuoa_test(fun, x0)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.solver = 'newuoa';

    [x, fval] = pdfo(problem);

end

function [x, fval] = uobyqa_test(fun, x0)

    problem = struct();
    problem.objective = fun;
    problem.x0 = x0;
    problem.solver = 'uobyqa';

    [x, fval] = pdfo(problem);

end