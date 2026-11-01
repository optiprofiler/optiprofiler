function options = validatePlotOptions(options)
%VALIDATEPLOTOPTIONS Native graphics property maps, not numerical options.
% Names are checked against the native classes. Values of every group except
% 'figure' are also checked before solving, by setting them on an unparented
% object of the same native class: no figure, window, renderer or
% process-global default is created or changed, and the installed MATLAB
% graphics API keeps its own validation. A figure object would open a
% window, so only the figure values whose rules need no figure are checked
% here (Position, OuterPosition, InnerPosition, PaperPosition, PaperSize,
% Color and Name), rejecting only values MATLAB rejects. Other figure values
% and handle-valued entries are checked when applied at export, where an
% error closes the figures of the call. Where graphics objects cannot be
% constructed (for example without a JVM), the values of the other groups are
% not checked here, while these figure values, which need no object, still
% are; portable output does not apply them.
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
        if strcmp(group, 'figure')
            checkFigureValues(options.(group));
        elseif ~isempty(names)
            checkValues(group, classes{index}, options.(group));
        end
    end
end

function checkFigureValues(values)
    names = fieldnames(values);
    for k = 1:numel(names)
        value = values.(names{k});
        % These properties also take the exact keywords of a settable default.
        if ismember(lower(names{k}), {'position', 'paperposition', 'papersize', 'color'}) && ...
                isScalarText(value) && any(strcmp(value, {'default', 'factory', 'remove'}))
            continue
        end
        switch lower(names{k})
            case 'position'
                valid = isRectangle(value) && all(real(value(3:4)) >= 0);
                rule = 'a vector of 4 finite numbers with nonnegative width and height';
            case {'outerposition', 'innerposition', 'paperposition'}
                valid = isRectangle(value);
                rule = 'a vector of 4 finite numbers';
            case 'papersize'
                valid = isnumeric(value) && numel(value) == 2 && ~any(real(value(:)) < 0);
                rule = 'a vector of 2 nonnegative numbers';
            case 'color'
                valid = isColor(value);
                rule = 'a color name, a hexadecimal color code, an RGB triplet or ''none''';
            case 'name'
                valid = ischar(value) || (isstring(value) && numel(value) <= 1 && ~any(ismissing(value)));
                rule = 'a character array or a string scalar';
            otherwise
                continue
        end
        if ~valid
            error('OptiProfiler:PlotOptions:InvalidProperty', ...
                'Invalid native plot_options.figure property: %s must be %s.', names{k}, rule);
        end
    end
end

function valid = isRectangle(value)
    valid = isnumeric(value) && numel(value) == 4 && all(isfinite(value(:)));
end

function valid = isColor(value)
    if (ischar(value) && isrow(value)) || (isstring(value) && isscalar(value))
        valid = strcmpi(value, 'none');
        if ~valid
            try
                validatecolor(value);
                valid = true;
            catch
            end
        end
    elseif isfloat(value) && isreal(value) && ismember(numel(value), [3, 4])
        valid = all(value(:) >= 0 & value(:) <= 1);
    else
        % Integer and other encodings are left to the figure.
        valid = ~(isnumeric(value) && ~ismember(numel(value), [3, 4])) && ~iscell(value);
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
