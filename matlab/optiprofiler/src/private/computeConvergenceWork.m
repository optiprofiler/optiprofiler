function [work_hist, work_out] = computeConvergenceWork(merit_histories, merit_outs, merit_inits, merit_mins, n_evals, tolerance, invalid_initial, valid_histories, valid_outs)
% Return first passing evaluation indices and work at passing solver outputs.
% Validity masks exclude raw NaN observations and padded history tails.
    n_problems = size(merit_histories, 1);
    n_solvers = size(merit_histories, 2);
    n_runs = size(merit_histories, 3);
    merit_init_inf_mask = isinf(merit_inits);
    work_hist = NaN(n_problems, n_solvers, n_runs);
    work_out = NaN(n_problems, n_solvers, n_runs);
    for i_problem = 1:n_problems
        for i_solver = 1:n_solvers
            for i_run = 1:n_runs
                if invalid_initial(i_problem, i_run)
                    continue;
                end
                if merit_init_inf_mask(i_problem, i_run)
                    % Preserve the Inf-initial convention only for
                    % actual evaluations without raw NaN values.
                    threshold = Inf;
                elseif isfinite(merit_mins(i_problem, i_run))
                    threshold = max(tolerance * merit_inits(i_problem, i_run) + (1 - tolerance) * merit_mins(i_problem, i_run), merit_mins(i_problem, i_run));
                else
                    % A -Inf minimum needs an exact -Inf hit; an
                    % undefined minimum cannot establish a target.
                    threshold = -Inf;
                end
                passing = valid_histories(i_problem, i_solver, i_run, :) & ...
                    (merit_histories(i_problem, i_solver, i_run, :) <= threshold);
                if any(passing(:))
                    work_hist(i_problem, i_solver, i_run) = find(passing, 1, 'first');
                end
                if valid_outs(i_problem, i_solver, i_run) && merit_outs(i_problem, i_solver, i_run) <= threshold
                    work_out(i_problem, i_solver, i_run) = n_evals(i_problem, i_solver, i_run);
                end
            end
        end
    end

end
