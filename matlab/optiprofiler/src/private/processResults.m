function [merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged] = processResults(results_plibs, profile_options)
%PROCESSRESULTS processes results from results_plibs.

    merit_histories_merged = [];
    merit_outs_merged = [];
    merit_inits_merged = [];
    n_evals_merged = [];
    problem_names_merged = [];
    problem_dims_merged = [];
    % Unify the length of merit_histories.
    results_plibs = unify_length(results_plibs, 'merit_histories');
    % Use results_plibs to merge the results from all the problem libraries.
    for i_plib = 1:size(results_plibs, 2)
        results_plib = results_plibs{i_plib};
        merit_histories_merged = cat(1, merit_histories_merged, results_plib.merit_histories);
        merit_outs_merged = cat(1, merit_outs_merged, results_plib.merit_outs);
        merit_inits_merged = cat(1, merit_inits_merged, results_plib.merit_inits);
        n_evals_merged = cat(1, n_evals_merged, results_plib.n_evals);
        problem_dims_merged = cat(1, problem_dims_merged, results_plib.problem_dims);
        problem_names_merged = [problem_names_merged, results_plib.problem_names];

        if isfield(results_plib, 'results_plib_plain')
            results_plib.merit_histories_plain = results_plib.results_plib_plain.merit_histories;
            results_plibs{i_plib} = results_plib;
        end
    end
    % Find the least merit value for each problem in each run.
    merit_mins_merged = min(min(merit_histories_merged, [], 4, 'omitnan'), [], 2, 'omitnan');
    for i_problem = 1:size(merit_histories_merged, 1)
        for i_run = 1:size(merit_histories_merged, 3)
            merit_mins_merged(i_problem, i_run) = min(merit_mins_merged(i_problem, i_run), merit_inits_merged(i_problem, i_run), 'omitnan');
        end
    end

    % If `results_plib_plain` exists and the `run_plain` field in `profile_options` is true. Then we need to redefine `merit_mins_merged`.
    if isfield(results_plibs{1}, 'results_plib_plain') && profile_options.(ProfileOptionKey.RUN_PLAIN.value)
        merit_histories_plain_merged = [];
        merit_inits_plain_merged = [];
        problem_names_plain_merged = [];
        % Unify the length of merit_histories_plain.
        results_plibs = unify_length(results_plibs, 'merit_histories_plain');
        for i_plib = 1:size(results_plibs, 2)
            results_plib = results_plibs{i_plib};
            merit_histories_plain_merged = cat(1, merit_histories_plain_merged, results_plib.merit_histories_plain);
            merit_inits_plain_merged = cat(1, merit_inits_plain_merged, results_plib.results_plib_plain.merit_inits);
            problem_names_plain_merged = [problem_names_plain_merged, results_plib.results_plib_plain.problem_names];
        end
        % Note that we will not consider merit_min for each run under the plain feature. Instead, we
        % only consider the minimum merit value among all runs under the plain feature.
        merit_mins_plain_merged = min(min(min(merit_histories_plain_merged, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
        for i_problem = 1:size(merit_histories_plain_merged, 1)
            merit_mins_plain_merged(i_problem) = min(merit_mins_plain_merged(i_problem), merit_inits_plain_merged(i_problem, i_run), 'omitnan');
        end
        for i_problem = 1:size(merit_histories_merged, 1)
            idx = find(strcmp(problem_names_merged{i_problem}, problem_names_plain_merged), 1);
            if isempty(idx)
                continue;
            end
            for i_run = 1:size(merit_histories_merged, 3)
                % Redefine the `merit_mins_merged` for the problems that are solved under the plain feature.
                % Note that min(x, NaN) = x.
                merit_mins_merged(i_problem, i_run) = min(merit_mins_merged(i_problem, i_run), merit_mins_plain_merged(idx), 'omitnan');
            end
            
        end
    end
end

function results_plibs = unify_length(results_plibs, name)
    % Unify the length of the specified field in the results_plibs which can be `merit_histories` or `merit_histories_plain`.

    max_length = 0;
    for i_plib = 1:size(results_plibs, 2)
        results_plib = results_plibs{i_plib};
        max_length = max(max_length, size(results_plib.(name), 4));
    end
    % Unify the length of the specified field in the results_plibs.
    for i_plib = 1:size(results_plibs, 2)
        results_plib = results_plibs{i_plib};
        field = results_plib.(name);
        field_unified = NaN(size(field, 1), size(field, 2), size(field, 3), max_length);
        field_unified(:, :, :, 1:size(field, 4)) = field;
        field_unified(:, :, :, size(field, 4) + 1:end) = repmat(field(:, :, :, end), [1, 1, 1, max_length - size(field, 4)]);
        results_plibs{i_plib}.(name) = field_unified;
    end
end