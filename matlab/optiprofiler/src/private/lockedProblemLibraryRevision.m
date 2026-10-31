function revision = lockedProblemLibraryRevision(name)
%LOCKEDPROBLEMLIBRARYREVISION returns the commit the integration lock names.
%   REVISION is the commit of NAME in matlab_problem_libraries.lock, or ''
%   when the lock has no such library or cannot be read. It is the reviewed
%   snapshot, not a verified identity of the installed checkout.

    revision = '';
    entry = lockedProblemLibraryEntry(name);
    if isstruct(entry) && isfield(entry, 'commit')
        revision = char(entry.commit);
    end
end
