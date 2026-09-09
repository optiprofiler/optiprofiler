function drawPerformanceDataProfiles(ax, x, y, solver_names, profile_options)
%DRAWPERFORMANCEDATAPROFILES draws performance profiles and data profiles.

    line_colors = profile_options.(ProfileOptionKey.LINE_COLORS.value);
    line_styles = profile_options.(ProfileOptionKey.LINE_STYLES.value);
    line_widths = profile_options.(ProfileOptionKey.LINE_WIDTHS.value);

    n_solvers = size(x, 2);
    [x_steps, mean_steps, lower_steps, upper_steps, n_runs] = prepareProfilePlotData(x, y, profile_options);

    hold(ax, 'on');

    for i_solver = 1:n_solvers
        x_stairs = x_steps{i_solver}; y_mean_stairs = mean_steps{i_solver};
        y_lower_stairs = lower_steps{i_solver}; y_upper_stairs = upper_steps{i_solver};

        % Get the color and the line style MATLAB will use for the next plot command in the axes 'ax'.
        color = line_colors(mod(i_solver - 1, size(line_colors, 1)) + 1, :);
        line_style = line_styles{mod(i_solver - 1, size(line_styles, 2)) + 1};
        line_width = line_widths(mod(i_solver - 1, length(line_widths)) + 1);

        % We first plot the shaded area and then the mean line to ensure that the shaded area is behind the mean line in case the installed MATLAB version does not support transparency.
        if n_runs > 1
            % The band is drawn as one patch of independent rectangular faces, one per interval of positive x-span and
            % positive height, instead of the single polygon [x; flip(x)], [lower; flip(upper)]: MATLAB tessellates a
            % patch face into triangles and, when the runs coincide on an interval and differ only at the duplicated x
            % of a jump, that self-touching polygon (zero geometric area at the jump) is painted as broad triangles
            % between its extreme corners. See optiprofiler_internal.profileBandFaces for the exact rule.
            [band_vertices, band_faces] = optiprofiler_internal.profileBandFaces(x_stairs, y_lower_stairs, y_upper_stairs);
            if ~isempty(band_faces)
                patch(ax, 'Vertices', band_vertices, 'Faces', band_faces, 'FaceColor', color, 'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        plot(ax, x_stairs, y_mean_stairs, line_style, 'Color', color, 'LineWidth', line_width, 'DisplayName', solver_names{i_solver});
    end

    set(ax, 'YLim', [0.0, 1.0]);
    set(ax, 'YTick', 0:0.2:1);
    yyaxis(ax, 'right');
    set(ax, 'YMinorTick','on');
    set(get(ax, 'YAxis'), 'MinorTickValues', [0.1 0.3 0.5 0.7 0.9]);
    set(ax, 'YTickLabel', []);
    yyaxis(ax, 'left');
    set(ax, 'YMinorTick','on');
    set(get(ax, 'YAxis'), 'MinorTickValues', [0.1 0.3 0.5 0.7 0.9]);
    linkprop([ax.XAxis; ax.YAxis],'color');
    linkprop([ax.YAxis(1), ax.YAxis(2)],{'Limits','TickValues'});
    box(ax,'on');
    hold(ax, 'off');
end
