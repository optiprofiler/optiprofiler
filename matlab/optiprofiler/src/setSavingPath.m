function [path_out, time_stamp] = setSavingPath(profile_options)
%SETSAVINGPATH Set the path to save the results and create the directory.

    % Paths to the results.
    % If 'benchmark_id' is not provided (by default '.'), we use time_stamp as the name of the directory.
    path_out = fullfile(profile_options.(ProfileOptionKey.SAVEPATH.value), 'out', profile_options.(ProfileOptionKey.BENCHMARK_ID.value));
    time_stamp = datestr(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH-mm-SSZ'), 'yyyy-mm-ddTHH-MM-SSZ');
    if strcmp(profile_options.(ProfileOptionKey.BENCHMARK_ID.value), '.')
        path_out = fullfile(path_out, time_stamp);
    end
    % Create the directory if it does not exist.
    if ~exist(path_out, 'dir')
        mkdir(path_out);
    end
end