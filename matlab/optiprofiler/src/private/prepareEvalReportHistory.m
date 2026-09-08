function panels = prepareEvalReportHistory(result, options)
%PREPAREEVALREPORTHISTORY Numeric history panels without graphics or callbacks.
%   Reuse the renderer's pure display/shift/aggregation pipeline. Do not
%   recompute merit: the collector can only use values the engine retained.
    panels = {};
    channels = {'objective', 'constraint', 'merit'};
    fields = {'fun_history', 'maxcv_history', 'merit_history'};
    initial_fields = {'fun_inits', 'maxcv_inits', 'merit_inits'};
    if strcmp(result.problem_type, 'u'), channels = channels(1); end
    for c = 1:numel(channels)
        if ~isfield(result, fields{c}) || ~isfield(result, initial_fields{c}), continue; end
        [display, note] = processHistYaxes(result.(fields{c}), result.(initial_fields{c}));
        shift = computeHistoryYShift(display, options);
        for is_cum = [false, true]
            [indices, center, lower, upper, n_runs] = prepareHistoryPlotData(display, is_cum, shift, result.n_eval, options);
            view = 'raw'; if is_cum, view = 'cummin'; end
            series = cell(1,numel(indices));
            for solver = 1:numel(indices)
                geometry = 'line'; if numel(indices{solver}) == 1, geometry = 'point'; end
                series{solver} = struct('solver_index', solver, ...
                    'evaluation_indices', {num2cell(indices{solver})}, ...
                    'x', {num2cell(indices{solver}/(result.problem_dim+1))}, ...
                    'mean', {num2cell(center{solver})}, 'lower', {num2cell(lower{solver})}, 'upper', {num2cell(upper{solver})}, ...
                    'geometry', geometry, 'band_visible', n_runs > 1 && numel(indices{solver}) > 1);
            end
            % One panel's scale is shared across all solvers, just like
            % drawFunMaxcvMeritHist (a varying nonzero curve enables log).
            scale = 'linear';
            if any(cellfun(@(v) any(v) && any(diff(v)), center)), scale = 'log'; end
            panels{end+1} = struct('kind', 'history', 'channel', channels{c}, 'mode', view, ...
                'status', 'numeric_prepared', 'fidelity', 'exact_rendered_data', 'series', {series}, ...
                'x_transform', 'evaluation_index/(dimension+1)', 'y_scale', scale, ...
                'y_shift', shift, 'display_limit', 1e100, 'display_note', note, ...
                'aggregation', options.hist_aggregation, 'errorbar_type', options.errorbar_type, ...
                'n_runs', n_runs, 'std_ddof', double(n_runs > 1));
        end
    end
end
