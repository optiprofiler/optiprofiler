function compareOneProblem()

    clc;
    solvers = {@bds_test, @fminsearch_test};
    labels = {'bds', 'fminsearch'};
    problem_name = 'arwhead';
    custom_problem_loader = @OPM_loader;
    feature_name = 'plain';
    feature_options = struct();
    testOneProblem(solvers, labels, problem_name, custom_problem_loader, feature_name, feature_options)

end


function x = bds_test(fun, x0)

    x = bds(fun, x0);

end

function x = fminsearch_test(fun, x0)

    x = fminsearch(fun, x0);

end

% function x = fminunc_test(fun, x0)
% 
%     % options = optimoptions("fminunc", "Algorithm", "quasi-newton", "HessUpdate", ...
%     %     fminunc_type, "MaxFunctionEvaluations", 1e20, "MaxIterations", 1e20, ...
%     %     "ObjectiveLimit", -Inf, "StepTolerance", eps, "OptimalityTolerance", eps, ...
%     %     'SpecifyObjectiveGradient', true);
%     x = fminunc(fun, x0);
% 
% end