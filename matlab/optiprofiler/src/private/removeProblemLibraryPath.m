function removeProblemLibraryPath(root)
%REMOVEPROBLEMLIBRARYPATH removes one registered root without deleting data.

    root = char(root);
    path_entries = strsplit(path, pathsep);
    matching = path_entries(strcmp(path_entries, root));
    if ~isempty(matching)
        original_warning_state = warning;
        warning('off', 'MATLAB:rmpath:DirNotFound');
        cleanup = onCleanup(@() warning(original_warning_state));
        rmpath(matching{:});
        clear cleanup
    end

    pathdef_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF');
    if isempty(pathdef_file)
        pathdef_file = fullfile(matlabroot, 'toolbox', 'local', 'pathdef.m');
    end
    original_warning_state = warning;
    warning('off', 'MATLAB:SavePath:PathNotSaved');
    cleanup = onCleanup(@() warning(original_warning_state));
    path_status = 1;
    path_message = '';
    try
        path_status = savepath(pathdef_file);
        if path_status ~= 0
            path_message = 'savepath returned a nonzero status.';
        end
    catch ME
        path_message = ME.message;
    end
    clear cleanup

    startup_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP');
    if isempty(startup_file) && ~isempty(userpath)
        startup_file = fullfile(userpath, 'startup.m');
    end
    try
        removeExactStartupEntry(startup_file, root);
    catch ME
        warning("MATLAB:unregisterProblemLibrary:startupCleanup", ...
            "The registration was removed, but startup path cleanup failed for '%s': %s", ...
            root, ME.message);
    end

    if path_status ~= 0
        warning("MATLAB:unregisterProblemLibrary:pathNotSaved", ...
            "The registration was removed, but MATLAB could not persist path removal for '%s': %s " + ...
            "Verify pathdef.m or startup.m before the next MATLAB session.", ...
            root, path_message);
    end
end


function removeExactStartupEntry(startup_file, root)
    if isempty(startup_file) || exist(startup_file, 'file') ~= 2
        return
    end

    contents = fileread(startup_file);
    lines = regexp(contents, '\r\n|\n|\r', 'split');
    expected = normalizeWhitespace(sprintf("addpath('%s'); %% optiprofiler", root));
    keep = true(size(lines));
    for i_line = 1:numel(lines)
        keep(i_line) = ~strcmp(normalizeWhitespace(lines{i_line}), expected);
    end
    if all(keep)
        return
    end

    updated = strjoin(lines(keep), newline);
    startup_dir = fileparts(startup_file);
    if isempty(startup_dir)
        startup_dir = pwd;
    end
    temporary_file = [tempname(startup_dir), '.m'];
    cleanup = onCleanup(@() deleteIfExists(temporary_file));
    file_id = fopen(temporary_file, 'w');
    if file_id == -1
        error('Cannot create temporary startup file %s.', temporary_file);
    end
    file_cleanup = onCleanup(@() fclose(file_id));
    fprintf(file_id, '%s', updated);
    clear file_cleanup
    [status, message] = movefile(temporary_file, startup_file, 'f');
    if ~status
        error('Cannot update startup file %s: %s', startup_file, message);
    end
    clear cleanup
end


function value = normalizeWhitespace(value)
    value = strtrim(regexprep(char(value), '\s+', ' '));
end


function deleteIfExists(filename)
    if exist(filename, 'file') == 2
        delete(filename);
    end
end
