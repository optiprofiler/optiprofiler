function entry = lockedProblemLibraryEntry(name)
%LOCKEDPROBLEMLIBRARYENTRY returns the entry of NAME in the integration lock.
%   ENTRY is the struct of NAME in matlab_problem_libraries.lock, or [] when
%   the lock has no such library or cannot be read.

    entry = [];
    try
        [lock, ~] = readMatlabProblemLibraryLock();
    catch
        return
    end
    entries = lock.libraries;
    for i_entry = 1:numel(entries)
        if iscell(entries)
            candidate = entries{i_entry};
        else
            candidate = entries(i_entry);
        end
        if isfield(candidate, 'name') && strcmp(char(candidate.name), name)
            entry = candidate;
            return
        end
    end
end
