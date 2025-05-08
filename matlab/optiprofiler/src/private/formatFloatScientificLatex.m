function [raw, formatted] = formatFloatScientificLatex(x, digits)
    % Format a floating-point number as scientific notation in LaTeX.

    if x < 0
        [raw, formatted] = formatFloatScientificLatex(-x, digits);
        raw = ['-', raw];
        formatted = ['-', formatted];
        return;
    end
    if x == 0
        raw = '0';
        formatted = '0';
        return;
    end

    exponent = floor(log10(abs(x)));
    coefficient = x / 10^exponent;
    % Round the coefficient to the specified number of digits.
    coefficient = round(coefficient, digits, "significant");

    % If the coefficient is 10, increase the exponent by 1.
    if coefficient == 10
        coefficient = 1;
        exponent = exponent + 1;
    end
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