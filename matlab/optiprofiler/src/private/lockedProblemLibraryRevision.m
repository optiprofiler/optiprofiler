function revision = lockedProblemLibraryRevision(name)
%LOCKEDPROBLEMLIBRARYREVISION returns the commit the integration lock names.
%   REVISION is the commit of NAME in matlab_problem_libraries.lock, or ''
%   when the lock has no such library or cannot be read. It is the reviewed
%   snapshot, not a verified identity of the installed checkout.

    revision = '';
    try
        [lock, ~] = readMatlabProblemLibraryLock();
    catch
        return
    end
    entries = lock.libraries;
    for i_entry = 1:numel(entries)
        if iscell(entries)
            entry = entries{i_entry};
        else
            entry = entries(i_entry);
        end
        if isfield(entry, 'name') && strcmp(char(entry.name), name) && isfield(entry, 'commit')
            revision = char(entry.commit);
            return
        end
    end
end
