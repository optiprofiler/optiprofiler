function [position, width, height] = profileFigurePosition(width, height)
%PROFILEFIGUREPOSITION returns the standard MATLAB profile figure position.
%
%   OptiProfiler's Python profiles use Matplotlib's default 6.4-by-4.8 inch
%   figure size. MATLAB figures are created with pixel units, so we use the
%   same 4:3 base size, 640-by-480, as the stable MATLAB counterpart.

    base_width = 640;
    base_height = 480;

    if nargin < 1 || isempty(width)
        width = base_width;
    end
    if nargin < 2 || isempty(height)
        height = base_height;
    end

    default_units = get(0, 'DefaultFigureUnits');
    if strcmpi(default_units, 'pixels')
        default_position = get(0, 'DefaultFigurePosition');
    else
        default_position = get(0, 'FactoryFigurePosition');
    end
    position = [default_position(1:2), width, height];
end
