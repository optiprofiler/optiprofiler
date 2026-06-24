function escaped = escapeLatexText(text)
%ESCAPELATEXTEXT escapes literal text for MATLAB's LaTeX interpreter.

    text = char(text);
    escaped = '';
    for i = 1:numel(text)
        switch text(i)
            case '\'
                piece = '\textbackslash{}';
            case '_'
                piece = '\_';
            case '#'
                piece = '\#';
            case '%'
                piece = '\%';
            case '&'
                piece = '\&';
            case '$'
                piece = '\$';
            case '{'
                piece = '\{';
            case '}'
                piece = '\}';
            case '^'
                piece = '\^{}';
            case '~'
                piece = '\~{}';
            otherwise
                piece = text(i);
        end
        escaped = [escaped, piece]; %#ok<AGROW>
    end
end
