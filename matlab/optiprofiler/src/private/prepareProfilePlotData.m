function [x_steps, mean_steps, lower_steps, upper_steps, n_runs] = prepareProfilePlotData(x, y, profile_options)
%PREPAREPROFILEPLOTDATA Shared native-profile and EvalReport presentation.
%   Preserve the renderer's aggregation and MATLAB sample std convention;
%   scores continue to use the unchanged unexpanded curves, never this data.
    n_solvers = size(x, 2);
    n_runs = size(y, 3);
    y_mean = squeeze(mean(y, 3));
    switch profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value)
        case 'minmax'
            y_lower = squeeze(min(y, [], 3));
            y_upper = squeeze(max(y, [], 3));
        case 'meanstd'
            y_std = squeeze(std(y, [], 3));
            y_lower = max(y_mean - y_std, 0);
            y_upper = min(y_mean + y_std, 1);
        otherwise
            error("Unknown `errorbar_type`: %s", profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value));
    end
    x_steps = cell(1,n_solvers); mean_steps = x_steps; lower_steps = x_steps; upper_steps = x_steps;
    for solver = 1:n_solvers
        [x_steps{solver}, mean_steps{solver}] = stairs(x(:,solver), y_mean(:,solver));
        [~, lower_steps{solver}] = stairs(x(:,solver), y_lower(:,solver));
        [~, upper_steps{solver}] = stairs(x(:,solver), y_upper(:,solver));
    end
end
