function library = resolveProblemLibrary(name)
%RESOLVEPROBLEMLIBRARY resolves and validates a MATLAB problem library.
%
%   LIBRARY = RESOLVEPROBLEMLIBRARY(NAME) returns a structure containing
%   canonical metadata plus validated SELECT and LOAD function handles.

    if ~ischarstr(name) || ~isscalar(string(name))
        error("MATLAB:resolveProblemLibrary:nameNotValid", ...
            "The problem-library name must be a char or string scalar.");
    end
    name = lower(strtrim(char(name)));
    if isempty(regexp(name, '^[a-z][a-z0-9_]*$', 'once'))
        error("MATLAB:resolveProblemLibrary:nameNotValid", ...
            "Invalid problem-library name '%s'.", name);
    end

    definitions = getProblemLibraryDefinitions();
    matching = find(strcmp({definitions.name}, name), 1);
    if isempty(matching)
        error("MATLAB:resolveProblemLibrary:notRegistered", ...
            "Problem library '%s' is not registered. " + ...
            "Use registerProblemLibrary to add external or custom libraries.", name);
    end
    registration = definitions(matching);
    if isempty(registration.root)
        error("MATLAB:resolveProblemLibrary:notInstalled", ...
            "Problem library '%s' is known but not installed. " + ...
            "Run setup or registerProblemLibrary first.", name);
    end

    current_platform = matlabPlatformName();
    if ~ismember(current_platform, registration.platforms)
        error("MATLAB:resolveProblemLibrary:unsupportedPlatform", ...
            "Problem library '%s' does not support MATLAB on %s.", ...
            name, current_platform);
    end
    library = materializeProblemLibrary(registration);
end
