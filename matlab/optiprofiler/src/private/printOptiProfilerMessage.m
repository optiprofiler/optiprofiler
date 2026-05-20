function printOptiProfilerMessage(level, message)
%PRINTOPTIPROFILERMESSAGE prints a wrapped OptiProfiler log message.

    lines = formatOptiProfilerMessage(level, message);
    for i_line = 1:numel(lines)
        fprintf('%s\n', lines{i_line});
    end
end

function lines = formatOptiProfilerMessage(level, message)
%FORMATOPTIPROFILERMESSAGE formats a wrapped OptiProfiler log message.

    width = 100;
    prefix = sprintf('%s: ', upper(char(level)));
    message = strtrim(regexprep(char(message), '\s+', ' '));
    if isempty(message)
        lines = {};
        return;
    end

    body_width = max(20, width - numel(prefix));
    body_lines = wrapTextForLog(message, body_width);
    continuation_prefix = repmat(' ', 1, numel(prefix));
    lines = cell(size(body_lines));
    lines{1} = [prefix, body_lines{1}];
    for i_line = 2:numel(body_lines)
        lines{i_line} = [continuation_prefix, body_lines{i_line}];
    end
end

function lines = wrapTextForLog(message, width)
    lines = {};
    remaining = message;
    while numel(remaining) > width
        chunk = remaining(1:width);
        break_idx = find(isspace(chunk), 1, 'last');
        if isempty(break_idx) || break_idx < 1
            break_idx = width;
        end
        lines{end + 1} = strtrim(remaining(1:break_idx)); %#ok<AGROW>
        remaining = strtrim(remaining(break_idx + 1:end));
    end
    lines{end + 1} = remaining;
end
