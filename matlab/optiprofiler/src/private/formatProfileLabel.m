function label = formatProfileLabel(template, tolerance, interpreter)
%FORMATPROFILELABEL substitutes the tolerance into a profile label template.
%   Only %s is replaced by TOLERANCE. %% and any other % give a literal percent
%   sign, written \% when INTERPRETER is 'latex', which refuses a bare % inside
%   and outside math mode, and copied as % otherwise. \% is copied as written,
%   also before s. Templates written for the earlier sprintf formatting keep
%   their meaning: \\ gives \, so \\%% gives \%. Every other character is
%   literal, including a LaTeX command such as \alpha. A string array is
%   formatted elementwise.

    if isstring(template)
        label = template;
        for k = 1:numel(template)
            label(k) = string(substitute(char(template(k)), tolerance, interpreter));
        end
    else
        label = substitute(template, tolerance, interpreter);
    end
end


function label = substitute(text, tolerance, interpreter)
    percent = '%';
    if strcmpi(interpreter, 'latex')
        percent = '\%';
    end
    % Longest token first at a backslash: \% before \\, and \\%% (or \\% that
    % is not \\ before the placeholder) as the one token it was for sprintf.
    [pieces, tokens] = regexp(text, '\\%|\\\\%(?:%|(?!s))|\\\\|%s|%%|%', 'split', 'match');
    replacements = tokens;
    replacements(strcmp(tokens, '%s')) = {char(tolerance)};
    replacements(strcmp(tokens, '%%') | strcmp(tokens, '%')) = {percent};
    replacements(strcmp(tokens, '\\')) = {'\'};
    replacements(strncmp(tokens, '\\%', 3)) = {'\%'};
    parts = [pieces; [replacements, {''}]];
    label = [parts{:}];
end
