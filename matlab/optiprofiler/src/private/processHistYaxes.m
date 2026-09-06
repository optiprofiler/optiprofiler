function [processed, display_note] = processHistYaxes(value_histories, value_inits)
%PROCESSHISTYAXES Prepare a bounded display copy, never experimental data.
%   The history shape is (solver, run, evaluation). Finite display values
%   outside +/-1e100 are clipped BEFORE computing means, squared deviations
%   and axis margins. This leaves ample double-precision arithmetic headroom.
%   NaN/Inf are displayed above their run's finite range, including its finite
%   initial value; with no finite reference they are displayed at 1.
%   Callers must show DISPLAY_NOTE when nonempty. This limit is not used by
%   oracle evaluations, saved histories, merit functions or profile scores.

    display_limit = 1e100;
    processed = max(-display_limit, min(display_limit, value_histories));
    notes = {};
    n_clipped = nnz(isfinite(value_histories) & abs(value_histories) > display_limit);
    n_nonfinite = nnz(~isfinite(value_histories));
    if n_clipped
        notes{end + 1} = sprintf('Display clipped at +/-1e100: %d entries', n_clipped);
    end
    if n_nonfinite
        notes{end + 1} = sprintf('Nonfinite placeholders: %d entries', n_nonfinite);
    end
    n_empty_runs = 0;
    for i_run = 1:size(value_histories, 2)
        % A slice mask must never index the whole multi-run array: MATLAB
        % linear indexing would otherwise change a different run's values.
        slice_run = processed(:, i_run, :);
        mask_run = ~isfinite(value_histories(:, i_run, :));
        if ~any(mask_run, 'all')
            continue;
        end
        finite = slice_run(~mask_run);
        initial = value_inits(min(i_run, numel(value_inits)));
        if isfinite(initial)
            finite(end + 1) = max(-display_limit, min(display_limit, initial));
        end
        if isempty(finite)
            replacement = 1;
            n_empty_runs = n_empty_runs + 1;
        else
            low = min(finite);
            high = max(finite);
            % Clip first so the range and placeholder cannot overflow. Give
            % a constant run a distinct, visible nonfinite placeholder too.
            gap = max([0.5 * (high - low), 0.05 * abs(high), eps]);
            replacement = high + gap;
        end
        slice_run(mask_run) = replacement;
        processed(:, i_run, :) = slice_run;
    end
    if n_empty_runs
        notes{end + 1} = sprintf('No finite reference in %d run(s): shown at 1', n_empty_runs);
    end
    display_note = strjoin(notes, newline);
end
