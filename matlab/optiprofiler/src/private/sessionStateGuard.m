function guard = sessionStateGuard()
%SESSIONSTATEGUARD Restores the session state that benchmark changes.
%   GUARD = SESSIONSTATEGUARD() records the caller's diary state, the root
%   defaults that benchmark sets and the MATLAB path. When GUARD is destroyed
%   (the caller returns or raises), it closes the experiment diary and
%   restores the caller's diary, restores those root defaults (removing a
%   default the caller had not set), removes the problem-library roots and the
%   folders under a library root that were added to the path in between, and
%   restores the order of the caller's entries, which adding a root that is
%   already on the path changes. Other path entries added in between, for
%   example by a solver, are kept. The guard MUST remain a local variable of
%   the caller.

    state.diary = get(0, 'Diary');
    state.diary_file = get(0, 'DiaryFile');
    defaults = get(groot, 'default');
    state.default_names = {'defaultAxesFontSize', 'defaultAxesFontName'};
    state.default_set = isfield(defaults, state.default_names);
    state.default_values = cell(size(state.default_names));
    for k = find(state.default_set)
        state.default_values{k} = defaults.(state.default_names{k});
    end
    state.path = strsplit(path, pathsep);
    guard = onCleanup(@() restore(state));
end


function restore(state)
    % Never propagate errors out of cleanup.
    try
        if ~strcmp(get(0, 'Diary'), state.diary) || ~strcmp(get(0, 'DiaryFile'), state.diary_file)
            diary('off');
            set(0, 'DiaryFile', state.diary_file);
            if strcmp(state.diary, 'on')
                diary('on');
            end
        end
    catch
    end
    for k = 1:numel(state.default_names)
        try
            if state.default_set(k)
                set(groot, state.default_names{k}, state.default_values{k});
            else
                set(groot, state.default_names{k}, 'remove');
            end
        catch
        end
    end
    try
        current = strsplit(path, pathsep);
        added = ~ismember(current, state.path);
        if any(added)
            kept = current(~added | ~libraryFolder(current));
        else
            kept = current;
        end
        % Solvers' entries keep their places; the caller's entries fill theirs.
        kept(ismember(kept, state.path)) = state.path(ismember(state.path, kept));
        if ~isequal(kept, current)
            path(strjoin(kept, pathsep));
        end
    catch
    end
end


function inside = libraryFolder(folders)
    % Whether each folder is a problem-library root or lies under one.
    inside = false(size(folders));
    try
        definitions = getProblemLibraryDefinitions();
    catch
        return
    end
    roots = {definitions.root};
    roots = roots(~cellfun(@isempty, roots));
    for k = 1:numel(roots)
        inside = inside | strcmp(folders, roots{k}) | startsWith(folders, [roots{k}, filesep]);
    end
end
