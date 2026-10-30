function exportGraphicsBounded(target, file, varargin)
%EXPORTGRAPHICSBOUNDED calls exportgraphics without waiting indefinitely.
%   Same arguments as exportgraphics. See boundedFigureCall for the bound.
    boundedFigureCall(ancestor(target, 'figure'), @() exportgraphics(target, file, varargin{:}));
end
