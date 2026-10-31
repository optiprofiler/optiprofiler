function label = formatProfileLabel(template, tolerance)
%FORMATPROFILELABEL substitutes the tolerance into a profile label template.
%   Only %s is replaced by TOLERANCE. Templates written for the earlier
%   sprintf formatting keep their meaning: %% gives % and \\ gives \. Every
%   other character is literal, including a LaTeX command such as \alpha and
%   a single %. A string array or a cell of chars is formatted elementwise.

    if isstring(template)
        label = template;
        for k = 1:numel(template)
            label(k) = string(substitute(char(template(k)), tolerance));
        end
    elseif iscell(template)
        label = cellfun(@(text) substitute(char(text), tolerance), template, 'UniformOutput', false);
    else
        label = substitute(template, tolerance);
    end
end


function label = substitute(text, tolerance)
    [pieces, tokens] = regexp(text, '%s|%%|\\\\', 'split', 'match');
    replacements = tokens;
    replacements(strcmp(tokens, '%s')) = {char(tolerance)};
    replacements(strcmp(tokens, '%%')) = {'%'};
    replacements(strcmp(tokens, '\\')) = {'\'};
    parts = [pieces; [replacements, {''}]];
    label = [parts{:}];
end
