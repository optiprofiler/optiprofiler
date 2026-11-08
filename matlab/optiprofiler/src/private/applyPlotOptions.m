function applyPlotOptions(figure_handle, profile_options)
%APPLYPLOTOPTIONS User properties win over drawing/layout defaults at export.
% No root defaults, objective callbacks, histories, or score arrays are changed.
    if ~isfield(profile_options,'plot_options') || isempty(fieldnames(profile_options.plot_options))
        return
    end
    options = profile_options.plot_options;
    apply(figure_handle,'figure');
    apply(findall(figure_handle,'Type','axes'),'axes');
    apply(findall(figure_handle,'Type','line'),'line');
    apply(findall(figure_handle,'Type','text'),'text');
    apply(findall(figure_handle,'Type','legend'),'legend');
    objects = findall(figure_handle);
    for k = 1:numel(objects)
        object = objects(k);
        if isa(object,'matlab.graphics.axis.Axes') || isa(object,'matlab.graphics.layout.TiledChartLayout')
            apply(object.XLabel,'xlabel');
            apply(object.YLabel,'ylabel');
            apply(object.Title,'title');
        end
        if isa(object,'matlab.graphics.axis.Axes')
            apply(object.XAxis,'xaxis');
            apply(object.YAxis,'yaxis');
        end
    end
    function apply(handles, group)
        if isempty(handles) || ~isfield(options,group) || isempty(fieldnames(options.(group))), return; end
        try
            set(handles,options.(group));
        catch cause
            error('OptiProfiler:PlotOptions:InvalidProperty','Invalid native plot_options.%s property: %s',group,cause.message);
        end
    end
end
