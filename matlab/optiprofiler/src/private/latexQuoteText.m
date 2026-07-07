function value = latexQuoteText(text)
%LATEXQUOTETEXT wraps escaped literal text in LaTeX double quotes.
%
%   Straight double quotes are typeset as closing quotes by MATLAB's LaTeX
%   interpreter, so use TeX's explicit opening and closing quote tokens.

    value = ['``', text, char([39 39])];
end
