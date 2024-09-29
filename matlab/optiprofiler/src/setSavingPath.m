function path_out = setSavingPath(profile_options)
%SETSAVINGPATH Set the path to save the results and create the directory.

    % Paths to the results.
    % If 'benchmark_id' is not provided (by default '.'), we use timestamp as the name of the directory.
    path_out = fullfile(profile_options.(ProfileOptionKey.SAVEPATH.value), 'out', profile_options.(ProfileOptionKey.BENCHMARK_ID.value));
    if strcmp(profile_options.(ProfileOptionKey.BENCHMARK_ID.value), '.')
        timestamp = datestr(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH-mm-SSZ'), 'yyyy-mm-ddTHH-MM-SSZ');
        path_out = fullfile(path_out, timestamp);
    end
end