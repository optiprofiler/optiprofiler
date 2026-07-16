function library = materializeProblemLibrary(registration)
%MATERIALIZEPROBLEMLIBRARY validates functions and returns callable handles.

    registration = normalizeProblemLibraryRegistration(registration);
    addpath(registration.root);
    rehash

    select_handle = resolveFunctionInRoot( ...
        registration.select_function, registration.root, true);
    load_handle = resolveFunctionInRoot( ...
        registration.load_function, registration.root, true);
    collect_info_handle = resolveFunctionInRoot( ...
        registration.collect_info_function, registration.root, false);
    check_available_handle = resolveFunctionInRoot( ...
        registration.check_available_function, registration.root, false);

    library = registration;
    library.select = select_handle;
    library.load = load_handle;
    library.collect_info = collect_info_handle;
    library.check_available = check_available_handle;
end


function handle = resolveFunctionInRoot(function_name, root, required)
    handle = [];
    if isempty(function_name)
        return
    end
    function_path = which(function_name);
    if isempty(function_path)
        if required
            error("MATLAB:problemLibraryRegistry:functionNotFound", ...
                "Required function '%s' was not found under '%s'.", ...
                function_name, root);
        end
        return
    end
    [status, attributes] = fileattrib(function_path);
    if status
        function_path = attributes.Name;
    end
    root_prefix = [root, filesep];
    if ~strcmp(fileparts(function_path), root) && ...
            ~startsWith(function_path, root_prefix)
        error("MATLAB:problemLibraryRegistry:functionOutsideRoot", ...
            "Function '%s' resolves outside registered root '%s': %s", ...
            function_name, root, function_path);
    end
    handle = str2func(function_name);
end
