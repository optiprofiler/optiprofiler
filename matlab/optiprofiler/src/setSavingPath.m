function path_out = setSavingPath(profile_options)
%SETSAVINGPATH Set the path to save the results and create the directory.

    % Paths to the results.
    timestamp = datestr(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH-mm-SSZ'), 'yyyy-mm-ddTHH-MM-SSZ');
    path_out = fullfile(profile_options.(ProfileOptionKey.SAVEPATH.value), 'out', profile_options.(ProfileOptionKey.BENCHMARK_ID.value));
    if ~exist(path_out, 'dir')
        mkdir(path_out);
    end

    % If the path does not exist or it is an empty directory, create it and store the timestamp, otherwise, append the timestamp.
    contents = dir(path_out);
    is_empty = isempty(contents) || (length(contents) == 2 && all(cellfun(@(x) ismember(x, {'.', '..'}), {contents.name})));

    if is_empty
        % Try to save the timestamp for later reference.
        try
            fid = fopen(fullfile(path_out, 'timestamp.txt'), 'w');
            fprintf(fid, timestamp);
            fclose(fid);
        catch
            fprintf("WARNING: The timestamp could not be saved.\n");
        end
    else
        % Check whether all the subdirectories are named by timestamps or 'time-unknown', 'time-unknown-2', 'time-unknown-3', ...
        timestamp_pattern = '^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}Z$';
        unknown_pattern = '^time-unknown(-\d+)?$';
        fullpattern = ['(' timestamp_pattern ')|(' unknown_pattern ')'];
        filtered_contents = contents(~startsWith({contents.name}, '.'));
        is_all_timestamps = all(arrayfun(@(x) x.isdir && ~isempty(regexp(x.name, fullpattern, 'once')), filtered_contents));
        if ~is_all_timestamps
            % Initialize the name of the directory by 'time-unknown'. If there are multiple directories named by 'time-unknown' or 'time-unknown-2', 'time-unknown-3', ..., use the next number.
            unknown_dirs = filtered_contents(contains({filtered_contents.name}, 'time-unknown'));
            if isempty(unknown_dirs)
                old_timestamp = 'time-unknown';
            else
                last_unknown_dir = unknown_dirs(end);
                last_unknown_dir_name = last_unknown_dir.name;
                if strcmp(last_unknown_dir_name, 'time-unknown')
                    old_timestamp = 'time-unknown-2';
                else
                    num_str = regexp(last_unknown_dir_name, '\d+$', 'match');
                    if isempty(num_str)
                        old_timestamp = 'time-unknown-2';
                    else
                        last_unknown_number = str2double(num_str{1});
                        old_timestamp = ['time-unknown-' num2str(last_unknown_number + 1)];
                    end
                end
            end
            % Try to find 'timestamp.txt' and read the timestamp.
            for i = 1:length(filtered_contents)
                if strcmp(filtered_contents(i).name, 'timestamp.txt')
                    try
                        fid = fopen(fullfile(path_out, 'timestamp.txt'), 'r');
                        old_timestamp_new = fscanf(fid, '%s');
                        fclose(fid);
                    catch
                    end
                    break;
                end
            end
            try
                mkdir(fullfile(path_out, old_timestamp_new));
                old_timestamp = old_timestamp_new;
            catch
                mkdir(fullfile(path_out, old_timestamp));
            end
            % Move all the filtered contents that does not match the fullpattern to the new directory.
            for i_item = 1:length(filtered_contents)
                item = filtered_contents(i_item);
                if isempty(regexp(item.name, fullpattern, 'once'))
                    movefile(fullfile(path_out, item.name), fullfile(path_out, old_timestamp, item.name));
                end
            end
        end
        path_out = fullfile(path_out, timestamp);
    end
    
    if ~exist(path_out, 'dir')
        mkdir(path_out);
    end
end