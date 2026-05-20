function short_message = shortenMessageForLog(message)
%SHORTENMESSAGEFORLOG shortens an error message for log output.

    max_length = 180;
    short_message = strtrim(regexprep(char(message), '\s+', ' '));
    if isempty(short_message)
        short_message = 'Unknown error.';
    elseif numel(short_message) > max_length
        short_message = [short_message(1:max_length-3), '...'];
    end
end
