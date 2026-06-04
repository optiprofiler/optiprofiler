function setProfileBoxAspect(ax)
%SETPROFILEBOXASPECT keeps profile plot boxes consistent across figures.

    defaultFigurePosition = get(0, 'DefaultFigurePosition');
    default_width = defaultFigurePosition(3);
    default_height = defaultFigurePosition(4);
    if default_width > eps
        pbaspect(ax, [default_width, default_height, 1]);
    end
end
