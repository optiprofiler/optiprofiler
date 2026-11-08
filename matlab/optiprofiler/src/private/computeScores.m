function profile_scores = computeScores(curves, profile_options)
    % Compute the scores of the solvers for all the profiles.

    n_tols = size(curves, 2);
    n_solvers = size(curves{1}.hist.perf, 1);

    if n_solvers == 2
        profile_scores = ones(n_solvers, n_tols, 2, 3);
    else
        profile_scores = ones(n_solvers, n_tols, 2, 2);
    end

    for i_tol = 1:n_tols
        for i_solver = 1:n_solvers
            curve_hist_perf = curves{i_tol}.hist.perf{i_solver, end};
            curve_hist_data = curves{i_tol}.hist.data{i_solver, end};
            curve_out_perf = curves{i_tol}.out.perf{i_solver, end};
            curve_out_data = curves{i_tol}.out.data{i_solver, end};
            profile_scores(i_solver, i_tol, 1, 1) = integrate(curve_hist_perf, 'perf', profile_options);
            profile_scores(i_solver, i_tol, 1, 2) = integrate(curve_hist_data, 'data', profile_options);
            profile_scores(i_solver, i_tol, 2, 1) = integrate(curve_out_perf, 'perf', profile_options);
            profile_scores(i_solver, i_tol, 2, 2) = integrate(curve_out_data, 'data', profile_options);
            if n_solvers == 2
                curve_hist_log_ratio = curves{i_tol}.hist.log_ratio{i_solver};
                curve_out_log_ratio = curves{i_tol}.out.log_ratio{i_solver};
                profile_scores(i_solver, i_tol, 1, 3) = integrate(curve_hist_log_ratio, 'log_ratio', profile_options);
                profile_scores(i_solver, i_tol, 2, 3) = integrate(curve_out_log_ratio, 'log_ratio', profile_options);
            end
        end
    end

    % Normalize the profile_scores by dividing the maximum profile_scores with the same tolerance and types.
    if profile_options.(ProfileOptionKey.NORMALIZED_SCORES.value)
        max_scores = max(profile_scores, [], 1);
        max_scores(max_scores == 0) = 1;
        profile_scores = profile_scores ./ max_scores;
    end
end
