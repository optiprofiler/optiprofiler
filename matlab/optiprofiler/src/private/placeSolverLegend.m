function placeSolverLegend(ax, n_solvers, default_location)
%PLACESOLVERLEGEND places the solver legend without hiding plotted curves.

    setProfileBoxAspect(ax);

    compact_threshold = 10;
    min_rows_per_column = 10;

    if n_solvers > compact_threshold
        n_cols = max(1, floor(n_solvers / min_rows_per_column));
        try
            lgd = legend(ax, 'Location', 'eastoutside', 'NumColumns', n_cols);
        catch
            lgd = legend(ax, 'Location', 'eastoutside');
        end
        if n_solvers <= 20
            lgd.FontSize = 8;
        elseif n_solvers <= 30
            lgd.FontSize = 7;
        else
            lgd.FontSize = 6;
        end
    else
        legend(ax, 'Location', chooseLegendLocation(ax, default_location));
    end
end


function location = chooseLegendLocation(ax, default_location)
    targets = struct(...
        'northeast', [0.85, 0.85], ...
        'northwest', [0.15, 0.85], ...
        'southeast', [0.85, 0.15], ...
        'southwest', [0.15, 0.15]);
    locations = fieldnames(targets);

    x_all = [];
    y_all = [];
    children = findobj(ax, 'Type', 'Line');
    for i_child = 1:numel(children)
        x_norm = normalizeAxisValues(children(i_child).XData, ax.XLim, ax.XScale);
        y_norm = normalizeAxisValues(children(i_child).YData, ax.YLim, ax.YScale);
        n = min(numel(x_norm), numel(y_norm));
        if n > 0
            x_all = [x_all, x_norm(1:n)]; %#ok<AGROW>
            y_all = [y_all, y_norm(1:n)]; %#ok<AGROW>
        end
    end

    if isempty(x_all)
        location = default_location;
        return;
    end

    best_score = Inf;
    location = default_location;
    for i_loc = 1:numel(locations)
        loc = locations{i_loc};
        target = targets.(loc);
        distance2 = (x_all - target(1)).^2 + (y_all - target(2)).^2;
        score = sum(exp(-distance2 / 0.06));
        if score < best_score || (abs(score - best_score) < eps && strcmp(loc, default_location))
            best_score = score;
            location = loc;
        end
    end
end


function normalized = normalizeAxisValues(values, limits, axis_scale)
    values = double(values);
    limits = double(limits);
    finite_mask = isfinite(values);
    if strcmp(axis_scale, 'log')
        finite_mask = finite_mask & values > 0 & all(limits > 0);
        if ~any(finite_mask)
            normalized = [];
            return;
        end
        values = log10(values(finite_mask));
        limits = log10(limits);
    else
        values = values(finite_mask);
    end
    span = limits(2) - limits(1);
    if ~isfinite(span) || abs(span) < eps
        normalized = [];
        return;
    end
    normalized = (values - limits(1)) / span;
    normalized = min(max(normalized(isfinite(normalized)), 0), 1);
    normalized = normalized(:)';
end
