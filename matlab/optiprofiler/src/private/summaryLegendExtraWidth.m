function extra_width = summaryLegendExtraWidth(n_solvers, default_width, solver_names)
%SUMMARYLEGENDEXTRAWIDTH extra per-profile summary width for outside legends.

    compact_threshold = 10;
    if n_solvers <= compact_threshold
        extra_width = 0;
        return;
    end
    if nargin < 3 || isempty(solver_names)
        max_label_length = 12;
    else
        max_label_length = max(cellfun(@(s) strlength(strrep(char(s), '\_', '_')), solver_names));
    end

    n_cols = legendColumnCount(n_solvers);
    per_column_fraction = min(0.45, max(0.22, 0.08 + 0.015 * double(max_label_length)));
    extra_width = default_width * (0.05 + per_column_fraction * n_cols);
end


function n_cols = legendColumnCount(n_solvers)
    min_rows_per_column = 10;
    n_cols = max(1, floor(n_solvers / min_rows_per_column));
end
