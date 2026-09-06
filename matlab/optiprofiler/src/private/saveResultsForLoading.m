function saveResultsForLoading(results_plibs, path_log, time_stamp, path_readme)
%SAVERESULTSFORLOADING Preserve numerical data before rendering can fail.
%   Version 7 may omit variables above 2 GB after only a warning. Stage a
%   version-7.3 payload, check its variable directory, then publish the marker.
%   A requested save failure is an error, including in silent mode.
    payload = fullfile(path_log, 'data_for_loading.mat');
    staging = fullfile(path_log, '.data_for_loading.mat');
    try
        save(staging, 'results_plibs', '-v7.3');
        if isempty(whos('-file', staging, 'results_plibs'))
            error('OptiProfiler:SavedResultsMissing', 'The saved file does not contain results_plibs.');
        end
    catch cause
        % An incomplete payload is not a recovery checkpoint. A COMPLETE
        % payload is retained if the subsequent atomic publication fails.
        if isfile(staging), delete(staging); end
        rethrow(cause);
    end
    atomicReplaceFile(staging, payload);
    marker_name = ['time_stamp_', time_stamp, '.txt'];
    marker = fullfile(path_log, marker_name);
    marker_staging = fullfile(path_log, '.time_stamp.txt');
    fid = fopen(marker_staging, 'w');
    if fid < 0
        error('OptiProfiler:TimeStamp', 'Cannot create the experiment time stamp file.');
    end
    cleanup = onCleanup(@() closeIfOpen(fid));
    written = fprintf(fid, '%s', time_stamp);
    closed = fclose(fid);
    clear cleanup;
    if written ~= numel(time_stamp) || closed ~= 0
        error('OptiProfiler:TimeStamp', 'Cannot finish writing the experiment time stamp file.');
    end
    atomicReplaceFile(marker_staging, marker);
    addToReadme(path_readme, 'data_for_loading.mat', 'File, storing the data of the current experiment for future loading.');
    addToReadme(path_readme, marker_name, 'File, recording the time stamp of the saved experiment.');
end

function closeIfOpen(fid)
    if ~isempty(fopen(fid)), fclose(fid); end
end
