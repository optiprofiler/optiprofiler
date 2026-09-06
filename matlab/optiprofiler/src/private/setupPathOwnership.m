function result = setupPathOwnership(action, context, owner, requested)
%SETUPPATHOWNERSHIP private records for paths actually introduced by setup.
% The public problem-library registry intentionally remains unchanged.
% POSIX startup/pathdef/sidecar persistence preserves existing modes; new files
% use 0600. Sensitive
% bytes are staged inside an atomically created 0700 directory before writing.
% Windows uses its existing ACL behavior, without a POSIX permission claim.
    filename = ownership_file();
    records = read_records(filename);
    result = false;
    if strcmp(action, 'add')
        context = char(context);
        owner = char(owner);
        requested = cellfun(@char, requested, 'UniformOutput', false);
        if ~absolute_path(context) || ~absolute_path(owner) || ...
                ~all(cellfun(@absolute_path, requested))
            error('OptiProfiler:PathOwnership', 'Setup ownership requires absolute paths.');
        end
        matching = find(strcmp({records.context}, context) & strcmp({records.owner}, owner));
        if isempty(matching)
            record = struct('context', context, 'owner', owner, ...
                'paths', {{}}, 'owned_paths', {{}}, 'startup_additions', {{}}, ...
                'pathdef_file', target_file('PATHDEF'), 'startup_file', target_file('STARTUP'), ...
                'pathdef_saved', false);
            records(end + 1) = record;
            matching = numel(records);
        elseif numel(matching) ~= 1
            error('OptiProfiler:PathOwnership', 'Duplicate setup ownership records.');
        end
        record = records(matching);
        if ~strcmp(record.pathdef_file, target_file('PATHDEF')) || ...
                ~strcmp(record.startup_file, target_file('STARTUP'))
            error('OptiProfiler:PathOwnership', ...
                'Persistence targets changed. Remove the previous setup ownership first.');
        end
        current = strsplit(path, pathsep);
        borrowed = current;
        if isfile(record.pathdef_file)
            borrowed = [borrowed, strsplit(read_persisted_path(record.pathdef_file), pathsep)];
        end
        for k = 1:numel(requested)
            entry = requested{k};
            if ~isfolder(entry)
                error('OptiProfiler:PathOwnership', 'Required path does not exist: %s', entry);
            end
            if ~ismember(entry, borrowed) && ~ismember(entry, record.owned_paths)
                record.owned_paths{end + 1} = entry;
            end
        end
        record.paths = unique([record.paths, requested], 'stable');
        records(matching) = record;
        write_records(filename, records);
        for k = 1:numel(requested)
            if ~ismember(requested{k}, current)
                addpath(requested{k});
            end
        end
        if isempty(record.owned_paths)
            result = true(numel(requested), 1);
            return;
        end
        saved = try_save_path(record.pathdef_file);
        record.pathdef_saved = record.pathdef_saved || saved;
        if ~saved
            if isempty(record.startup_file)
                error('OptiProfiler:PathPersistence', 'No writable pathdef or explicit startup fallback.');
            end
            for k = 1:numel(record.owned_paths)
                entry = record.owned_paths{k};
                quoted = strrep(entry, '''', '''''');
                addition = sprintf('\naddpath(''%s''); %% optiprofiler setup-owned', quoted);
                existing = read_bytes(record.startup_file);
                if ~contains(char(existing), addition)
                    % A completed write is recorded only after publication.
                    % An identical unrecorded line is borrowed, never claimed.
                    append_bytes(record.startup_file, unicode2native(addition, 'UTF-8'));
                elseif ~ismember(addition, record.startup_additions)
                    warning('OptiProfiler:UnownedStartupEntry', ...
                        'An existing startup entry for %s has no completed ownership record; it will not be deleted.', entry);
                    continue;
                end
                if ~ismember(addition, record.startup_additions)
                    record.startup_additions{end + 1} = addition;
                    records(matching) = record;
                    write_records(filename, records);
                end
            end
        end
        records(matching) = record;
        write_records(filename, records);
        result = true(numel(requested), 1);
        return;
    elseif strcmp(action, 'remove-context')
        matching = find(strcmp({records.context}, char(context)));
    elseif strcmp(action, 'remove-owner')
        matching = find(strcmp({records.owner}, char(context)));
    else
        error('OptiProfiler:PathOwnership', 'Unknown path ownership operation.');
    end
    result = ~isempty(matching);
    for index = matching
        record = records(index);
        current = strsplit(path, pathsep);
        remove = record.owned_paths(ismember(record.owned_paths, current));
        if ~isempty(remove)
            rmpath(remove{:});
        end
        if record.pathdef_saved && ~try_save_path(record.pathdef_file, record.owned_paths)
            error('OptiProfiler:PathPersistence', ...
                'Could not persist removal to %s; ownership retained for retry.', record.pathdef_file);
        end
        remove_additions(record.startup_file, record.startup_additions);
        % Keep a tombstone so a later unregister does not delete borrowed paths.
        records(index).owned_paths = {};
        records(index).startup_additions = {};
        records(index).pathdef_saved = false;
        write_records(filename, records);
    end
end

function filename = ownership_file()
    registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    if isempty(registry)
        registry = fullfile(prefdir, 'optiprofiler', 'problem_libraries.mat');
    end
    filename = [registry, '.setup-paths.mat'];
    if ~absolute_path(filename)
        error('OptiProfiler:PathOwnership', 'The registry/ownership location must be absolute.');
    end
end

function filename = target_file(kind)
    filename = getenv(['OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_', kind]);
    if isempty(filename)
        if strcmp(kind, 'PATHDEF')
            filename = fullfile(matlabroot, 'toolbox', 'local', 'pathdef.m');
        elseif ~isempty(userpath)
            filename = fullfile(userpath, 'startup.m');
        end
    end
    if ~isempty(filename) && ~absolute_path(filename)
        error('OptiProfiler:PathOwnership', 'Path persistence targets must be absolute.');
    end
end

function valid = absolute_path(value)
    valid = ischar(value) && isrow(value) && ~isempty(value) && ...
        ~any(ismember(value, [char(10), char(13)]));
    if ~valid
        return;
    end
    if ispc
        % A single leading separator is relative to the current drive.
        % Match setup.m: require a drive root or a UNC prefix on Windows.
        valid = ~isempty(regexp(value, '^[A-Za-z]:[\\/]', 'once')) || ...
            startsWith(value, '\\');
    else
        valid = startsWith(value, filesep);
    end
end

function records = read_records(filename)
    records = struct('context', {}, 'owner', {}, 'paths', {}, 'owned_paths', {}, ...
        'startup_additions', {}, 'pathdef_file', {}, 'startup_file', {}, 'pathdef_saved', {});
    if ~isfile(filename)
        return;
    end
    data = load(filename, 'schema_version', 'records');
    if ~isfield(data, 'schema_version') || ~isequal(data.schema_version, 1) || ...
            ~isfield(data, 'records') || ~isstruct(data.records) || ...
            ~all(isfield(data.records, fieldnames(records)))
        error('OptiProfiler:PathOwnership', 'Invalid setup ownership file: %s', filename);
    end
    records = data.records;
end

function write_records(filename, records)
    directory = fileparts(filename);
    if ~isfolder(directory)
        mkdir(directory);
    end
    [temporary, staging, mode] = private_stage(filename);
    schema_version = 1;
    try
        save(temporary, 'schema_version', 'records');
    catch cause
        discard_stage(staging);
        rethrow(cause);
    end
    publish_persistence(temporary, filename, staging, mode);
end

function saved = try_save_path(filename, removed)
    reject_symlink(filename);
    if nargin < 2
        removed = {};
    end
    state = warning;
    cleanup = onCleanup(@() warning(state)); %#ok<NASGU>
    warning('off', 'MATLAB:SavePath:PathNotSaved');
    if isfile(filename)
        % Preserve paths already persisted by the user even when they are
        % absent from this session. Remove only the recorded owned entries.
        original_path = path;
        path_cleanup = onCleanup(@() path(original_path)); %#ok<NASGU>
        persisted = read_persisted_path(filename);
        entries = unique([strsplit(persisted, pathsep), strsplit(original_path, pathsep)], 'stable');
        entries = entries(~ismember(entries, removed) & ~cellfun(@isempty, entries));
        path(strjoin(entries, pathsep));
    end
    directory = fileparts(filename);
    saved = false;
    if ~isfolder(directory)
        return;
    end
    if isfile(filename)
        [~, attributes] = fileattrib(filename);
        if ~attributes.UserWrite
            return;
        end
    end
    [temporary, staging, mode] = private_stage(filename);
    try
        status = savepath(temporary);
    catch cause
        discard_stage(staging);
        rethrow(cause);
    end
    if status ~= 0
        discard_stage(staging);
        return;
    end
    publish_persistence(temporary, filename, staging, mode);
    saved = true;
end

function persisted = read_persisted_path(filename)
    [directory, name, extension] = fileparts(filename);
    if ~strcmp(extension, '.m') || ~isvarname(name)
        error('OptiProfiler:PathPersistence', 'Cannot read persisted MATLAB path: %s', filename);
    end
    original_dir = pwd;
    cleanup = onCleanup(@() cd(original_dir)); %#ok<NASGU>
    cd(directory);
    clear(name);
    persisted = feval(name);
    if ~ischar(persisted) || ~isrow(persisted)
        error('OptiProfiler:PathPersistence', 'Invalid persisted MATLAB path in %s.', filename);
    end
end

function bytes = read_bytes(filename)
    bytes = uint8([]);
    if ~isfile(filename)
        return;
    end
    fid = fopen(filename, 'rb');
    if fid < 0
        error('OptiProfiler:PathPersistence', 'Cannot read %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    bytes = fread(fid, Inf, '*uint8').';
end

function append_bytes(filename, bytes)
    original = read_bytes(filename);
    publish_startup(filename, original, [original, bytes]);
end

function publish_startup(filename, original, updated)
    directory = fileparts(filename);
    if ~isfolder(directory)
        error('OptiProfiler:PathPersistence', 'Startup parent directory does not exist: %s.', directory);
    end
    if isfile(filename)
        [~, attributes] = fileattrib(filename);
        if ~attributes.UserWrite
            error('OptiProfiler:PathPersistence', 'Startup file is not writable: %s.', filename);
        end
    end
    [temporary, staging, mode] = private_stage(filename);
    try
        fid = fopen(temporary, 'wb');
        if fid < 0
            error('OptiProfiler:PathPersistence', 'Cannot create temporary startup file.');
        end
        file_cleanup = onCleanup(@() fclose(fid));
        if fwrite(fid, updated, 'uint8') ~= numel(updated)
            error('OptiProfiler:PathPersistence', 'Incomplete startup write.');
        end
        clear file_cleanup
    catch cause
        clear file_cleanup
        discard_stage(staging);
        rethrow(cause);
    end
    if ~isequal(read_bytes(filename), original)
        error('OptiProfiler:PathPersistence', ...
            'Startup changed during setup; original retained, candidate at %s.', temporary);
    end
    publish_persistence(temporary, filename, staging, mode);
end

function [temporary, staging, mode] = private_stage(filename)
    reject_symlink(filename);
    parent = fileparts(filename);
    mode = '';
    if isunix
        mode = '600';
        if isfile(filename)
            if ismac
                command = '/usr/bin/stat -L -f %Lp ';
            else
                command = '/usr/bin/stat -L -c %a ';
            end
            [status, output] = system([command, quotePosix(filename), ' 2>&1']);
            mode = strtrim(output);
            if status ~= 0 || isempty(regexp(mode, '^[0-7]{1,4}$', 'once'))
                error('OptiProfiler:PathPersistence', 'Cannot read POSIX mode for %s.', filename);
            end
        end
        % mktemp creates 0700 atomically, so no sensitive byte is ever written
        % in a directory with the caller's potentially wider default umask.
        template = fullfile(parent, '.optiprofiler-setup-XXXXXXXX');
        [status, output] = system(['/usr/bin/mktemp -d ', quotePosix(template), ' 2>&1']);
        if status ~= 0
            error('OptiProfiler:PathPersistence', 'Cannot create private staging directory: %s', strtrim(output));
        end
        staging = strtrim(output);
        [ok, attributes] = fileattrib(staging);
        if ~ok || ~attributes.UserRead || ~attributes.UserWrite || ~attributes.UserExecute || ...
                attributes.GroupRead || attributes.GroupWrite || attributes.GroupExecute || ...
                attributes.OtherRead || attributes.OtherWrite || attributes.OtherExecute
            error('OptiProfiler:PathPersistence', 'Staging directory is not private: %s.', staging);
        end
    else
        staging = tempname(parent);
        [ok, message] = mkdir(staging);
        if ~ok
            error('OptiProfiler:PathPersistence', 'Cannot create staging directory: %s', message);
        end
    end
    [~, ~, extension] = fileparts(filename);
    temporary = fullfile(staging, ['candidate', extension]);
end

function reject_symlink(filename)
    % Rename would replace the link itself, silently disconnecting a dotfile.
    % Check only the final component (including dangling links), not parents.
    if isunix
        status = system(['/bin/test -L ', quotePosix(filename)]);
        if status == 0
            error('OptiProfiler:SymlinkPersistence', ...
                'Refusing to replace symbolic-link persistence target %s. Choose an explicit regular-file target; the link was preserved.', filename);
        elseif status ~= 1
            error('OptiProfiler:PathPersistence', 'Cannot inspect persistence target %s.', filename);
        end
    end
end

function publish_persistence(temporary, filename, staging, mode)
    try
        if isunix
            [status, output] = system(['/bin/chmod ', mode, ' ', quotePosix(temporary), ' 2>&1']);
            if status ~= 0
                error('OptiProfiler:PathPersistence', 'Cannot set candidate POSIX mode: %s', strtrim(output));
            end
        end
        atomicReplaceFile(temporary, filename);
    catch cause
        % Never delete a completed candidate after a publication failure.
        failure = MException('OptiProfiler:PathPersistence', ...
            'Persistence failed for %s; complete candidate retained inside private staging at %s.', filename, temporary);
        throw(addCause(failure, cause));
    end
    discard_stage(staging);
end

function discard_stage(staging)
    if isfolder(staging)
        rmdir(staging, 's');
    end
end

function remove_additions(filename, additions)
    if isempty(filename) || ~isfile(filename) || isempty(additions)
        return;
    end
    original = read_bytes(filename);
    updated = original;
    for k = 1:numel(additions)
        token = unicode2native(additions{k}, 'UTF-8');
        positions = strfind(char(updated), char(token));
        if numel(positions) > 1
            error('OptiProfiler:PathOwnership', ...
                'Ambiguous owned startup entry in %s; no startup bytes changed.', filename);
        elseif numel(positions) == 1
            updated(positions:positions + numel(token) - 1) = [];
        end
    end
    if isequal(updated, original)
        return;
    end
    publish_startup(filename, original, updated);
end
