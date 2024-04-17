function [raw, formatted] = formatFloatScientificLatex(x)
    % Format a floating-point number as scientific notation in LaTeX.

    if x == 0
        raw = '0';
        formatted = '0';
        return;
    end

    exponent = floor(log10(abs(x)));
    coefficient = x / 10^exponent;
    formattedCoefficient = sprintf('%g', coefficient);
    if exponent == 0
        formatted = formattedCoefficient;
        raw = formatted;
    elseif coefficient == 1
        formatted = sprintf('10^{%d}', exponent);
        raw = sprintf('%g', x);
    else
        raw = sprintf('%se%d', formattedCoefficient, exponent);
        formatted = sprintf('%s \\times 10^{%d}', formattedCoefficient, exponent);
    end
    
end