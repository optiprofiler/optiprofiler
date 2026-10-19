function options = validatePlotOptions(options)
%VALIDATEPLOTOPTIONS Native graphics property maps, not numerical options.
% Values retain the installed MATLAB graphics API's own validation. Inspect
% names without creating graphics or changing process-global defaults.
    groups = {'figure','axes','xaxis','yaxis','xlabel','ylabel','title','text','line','legend'};
    classes = {'matlab.ui.Figure','matlab.graphics.axis.Axes', ...
        'matlab.graphics.axis.decorator.NumericRuler','matlab.graphics.axis.decorator.NumericRuler', ...
        'matlab.graphics.primitive.Text','matlab.graphics.primitive.Text', ...
        'matlab.graphics.primitive.Text','matlab.graphics.primitive.Text', ...
        'matlab.graphics.chart.primitive.Line','matlab.graphics.illustration.Legend'};
    if ~(isstruct(options) && isscalar(options))
        error('OptiProfiler:PlotOptions:NotStruct','plot_options must be a scalar struct of native property structs.');
    end
    fields = fieldnames(options);
    for k = 1:numel(fields)
        group = fields{k};
        index = find(strcmp(group,groups),1);
        if isempty(index)
            error('OptiProfiler:PlotOptions:UnknownGroup','Unknown plot_options group: %s.',group);
        end
        if ~(isstruct(options.(group)) && isscalar(options.(group)))
            error('OptiProfiler:PlotOptions:NotStruct','plot_options.%s must be a scalar native property struct.',group);
        end
        names = fieldnames(options.(group));
        supported = properties(classes{index});
        for n = 1:numel(names)
            if ~any(strcmpi(names{n},supported))
                error('OptiProfiler:PlotOptions:UnknownProperty','Unknown native property plot_options.%s.%s.',group,names{n});
            end
        end
    end
end
