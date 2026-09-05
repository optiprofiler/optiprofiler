function key = parseExperimentTimeStamp(time_stamp)
%PARSEEXPERIMENTTIMESTAMP Sort old timestamps and optional collision suffixes.
%   An invalid identifier returns [NaN, NaN] and is ignored during discovery.

    key = [NaN, NaN];
    parts = regexp(char(time_stamp), '^(\d{8}_\d{6})(?:_(\d{3,}))?$', 'tokens', 'once');
    if isempty(parts)
        return;
    end
    sequence = 0;
    if numel(parts) > 1 && ~isempty(parts{2})
        sequence = str2double(parts{2});
        if ~isfinite(sequence) || sequence == 0
            return;
        end
    end
    try
        time = datetime(parts{1}, 'InputFormat', 'yyyyMMdd_HHmmss');
        key = [datenum(time), sequence];
    catch
        % Malformed legacy markers are not loadable experiments.
    end
end
