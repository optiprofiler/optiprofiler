function [path_out, path_benchmark, time_stamp] = setSavingPath(profile_options)
%SETSAVINGPATH Set the path to save the results and create the directory.

    % Paths to the results.
    % If 'benchmark_id' is not provided (by default '.'), we use time_stamp as the name of the directory.
    path_out = fullfile(profile_options.(ProfileOptionKey.SAVEPATH.value), 'out');
    path_benchmark = fullfile(path_out, profile_options.(ProfileOptionKey.BENCHMARK_ID.value));
    time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
    time_stamp = char(time);
    time_stamp_folder = char(datetime(time, 'Format', 'yyyy-MM-dd''T''HH-mm-SS'));
    if strcmp(profile_options.(ProfileOptionKey.BENCHMARK_ID.value), '.')
        path_benchmark = fullfile(path_benchmark, time_stamp_folder);
    end
    % Create the directory if it does not exist.
    if ~exist(path_benchmark, 'dir')
        mkdir(path_benchmark);
    end
end