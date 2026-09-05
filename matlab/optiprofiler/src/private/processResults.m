function [merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged] = processResults(results_plibs, profile_options)
%PROCESSRESULTS processes results from results_plibs.

    merit_histories_merged = [];
    merit_outs_merged = [];
    merit_inits_merged = [];
    n_evals_merged = [];
    problem_names_merged = [];
    problem_dims_merged = [];
    % Old saved custom-merit arrays may predate raw-NaN validation.
    for i_plib = 1:numel(results_plibs)
        results_plibs{i_plib} = maskInvalidMerits(results_plibs{i_plib});
    end
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
    end
    % Find the least merit value for each problem in each run.
    merit_mins_merged = min(min(merit_histories_merged, [], 4, 'omitnan'), [], 2, 'omitnan');
    for i_problem = 1:size(merit_histories_merged, 1)
        for i_run = 1:size(merit_histories_merged, 3)
            merit_mins_merged(i_problem, i_run) = min(merit_mins_merged(i_problem, i_run), merit_inits_merged(i_problem, i_run), 'omitnan');
        end
    end

    % Match plain rows inside their original library container. MAT files
    % already preserve this nesting, including after load-time filtering.
    if isfield(profile_options, ProfileOptionKey.RUN_PLAIN.value) && profile_options.(ProfileOptionKey.RUN_PLAIN.value)
        first_problem = 0;
        for i_plib = 1:size(results_plibs, 2)
            results_plib = results_plibs{i_plib};
            names = results_plib.problem_names;
            offset = first_problem;
            first_problem = first_problem + numel(names);
            if ~isfield(results_plib, 'results_plib_plain') || isempty(results_plib.results_plib_plain)
                continue;
            end
            plain = results_plib.results_plib_plain;
            plain_names = plain.problem_names;
            if isempty(plain_names)
                continue;
            end
            if numel(unique(names)) ~= numel(names) || numel(unique(plain_names)) ~= numel(plain_names)
                error('MATLAB:processResults:duplicateProblemNames', ...
                    'The run_plain reference is ambiguous: duplicate problem names within one problem library cannot be matched uniquely.');
            end
            % One reference per problem, over all plain solvers/runs/evaluations.
            plain_mins = min(min(min(plain.merit_histories, [], 4, 'omitnan'), [], 3, 'omitnan'), [], 2, 'omitnan');
            for i_plain = 1:numel(plain_names)
                % Preserve the existing first-run convention for plain inits.
                plain_mins(i_plain) = min(plain_mins(i_plain), plain.merit_inits(i_plain), 'omitnan');
            end
            for i_problem = 1:numel(names)
                idx = find(strcmp(names{i_problem}, plain_names), 1);
                if isempty(idx)
                    continue;
                end
                row = offset + i_problem;
                merit_mins_merged(row, :) = min(merit_mins_merged(row, :), plain_mins(idx), 'omitnan');
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
