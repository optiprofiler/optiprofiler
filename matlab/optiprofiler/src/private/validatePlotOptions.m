function options = validatePlotOptions(options)
%VALIDATEPLOTOPTIONS Native graphics property maps, not numerical options.
% Names are checked against the native classes. Values of every group except
% 'figure' are also checked before solving, by setting them on an unparented
% object of the same native class: no figure, window, renderer or
% process-global default is created or changed, and the installed MATLAB
% graphics API keeps its own validation. Figure properties (a figure object
% would open a window) and handle-valued entries are checked when applied at
% export. Where graphics objects cannot be constructed (for example without a
% JVM), values are not checked here; portable output does not apply them.
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
        if ~strcmp(group, 'figure') && ~isempty(names)
            checkValues(group, classes{index}, options.(group));
        end
    end
end

function checkValues(group, class_name, values)
    names = fieldnames(values);
    % Object references (handles, parents, menus, layouts) could reparent or
    % alias the probe; leave them to the export-time application, which
    % reports the same error identifier.
    relational = {'Parent', 'Children', 'UIContextMenu', 'ContextMenu', 'Layout', 'CurrentPoint'};
    keep = cellfun(@(name) ~any(strcmpi(name, relational)) && ...
        ~(isobject(values.(name)) && ~isstring(values.(name))), names);
    values = rmfield(values, names(~keep));
    if isempty(fieldnames(values))
        return
    end
    try
        probe = feval(class_name);
    catch
        return
    end
    cleanup = onCleanup(@() delete(probe));
    try
        set(probe, values);
    catch cause
        error('OptiProfiler:PlotOptions:InvalidProperty','Invalid native plot_options.%s property: %s',group,cause.message);
    end
end
