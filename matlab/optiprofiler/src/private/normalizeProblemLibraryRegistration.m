function registration = normalizeProblemLibraryRegistration(registration)
%NORMALIZEPROBLEMLIBRARYREGISTRATION validates registry metadata.

    if ~isstruct(registration) || ~isscalar(registration)
        error("MATLAB:problemLibraryRegistry:registrationNotValid", ...
            "A problem-library registration must be a scalar structure.");
    end
    required = {'name', 'root', 'select_function', 'load_function'};
    missing = required(~isfield(registration, required));
    if ~isempty(missing)
        error("MATLAB:problemLibraryRegistry:registrationNotValid", ...
            "The problem-library registration is missing: %s.", strjoin(missing, ', '));
    end

    registration.name = normalizeName(registration.name, 'name');
    registration.root = normalizeString(registration.root, 'root', false);
    if exist(registration.root, 'dir') ~= 7
        error("MATLAB:problemLibraryRegistry:rootNotFound", ...
            "The problem-library root '%s' does not exist.", registration.root);
    end
    [status, attributes] = fileattrib(registration.root);
    if status
        registration.root = attributes.Name;
    end
    registration.select_function = normalizeName(registration.select_function, 'select_function');
    registration.load_function = normalizeName(registration.load_function, 'load_function');

    if ~isfield(registration, 'api_version')
        registration.api_version = 1;
    end
    if ~isintegerscalar(registration.api_version) || registration.api_version ~= 1
        error("MATLAB:problemLibraryRegistry:apiVersionNotSupported", ...
            "MATLAB problem libraries must use API version 1.");
    end
    registration.api_version = double(registration.api_version);

    if ~isfield(registration, 'role')
        registration.role = 'external';
    end
    registration.role = normalizeString(registration.role, 'role', false);

    if ~isfield(registration, 'collect_info_function') || isempty(registration.collect_info_function)
        registration.collect_info_function = '';
    else
        registration.collect_info_function = normalizeName( ...
            registration.collect_info_function, 'collect_info_function');
    end
    if ~isfield(registration, 'check_available_function') || isempty(registration.check_available_function)
        registration.check_available_function = '';
    else
        registration.check_available_function = normalizeName( ...
            registration.check_available_function, 'check_available_function');
    end

    if ~isfield(registration, 'platforms') || isempty(registration.platforms)
        registration.platforms = {'linux', 'macos', 'windows'};
    elseif ischarstr(registration.platforms)
        registration.platforms = {char(registration.platforms)};
    end
    if ~iscell(registration.platforms) || ...
            ~all(cellfun(@ischarstr, registration.platforms))
        error("MATLAB:problemLibraryRegistry:platformsNotValid", ...
            "The platforms field must be a nonempty cell array of names.");
    end
    registration.platforms = cellfun( ...
        @(platform) lower(char(platform)), registration.platforms, ...
        'UniformOutput', false);
    registration.platforms = unique(registration.platforms, 'stable');
    allowed_platforms = {'linux', 'macos', 'windows'};
    if isempty(registration.platforms) || ...
            ~all(ismember(registration.platforms, allowed_platforms))
        error("MATLAB:problemLibraryRegistry:platformsNotValid", ...
            "Platforms must be selected from linux, macos, and windows.");
    end

    ordered_fields = { ...
        'name', 'api_version', 'role', 'root', 'select_function', ...
        'load_function', 'collect_info_function', ...
        'check_available_function', 'platforms'};
    unknown = setdiff(fieldnames(registration), ordered_fields);
    if ~isempty(unknown)
        error("MATLAB:problemLibraryRegistry:registrationNotValid", ...
            "Unknown registration field(s): %s.", strjoin(unknown, ', '));
    end
    registration = orderfields(registration, ordered_fields);
end


function value = normalizeName(value, field_name)
    value = normalizeString(value, field_name, false);
    if isempty(regexp(value, '^[A-Za-z][A-Za-z0-9_]*$', 'once'))
        error("MATLAB:problemLibraryRegistry:registrationNotValid", ...
            "The field %s must be a MATLAB identifier.", field_name);
    end
    if strcmp(field_name, 'name')
        value = lower(value);
    end
end


function value = normalizeString(value, field_name, allow_empty)
    if ~ischarstr(value) || ~isscalar(string(value))
        error("MATLAB:problemLibraryRegistry:registrationNotValid", ...
            "The field %s must be a char or string scalar.", field_name);
    end
    value = strtrim(char(value));
    if ~allow_empty && isempty(value)
        error("MATLAB:problemLibraryRegistry:registrationNotValid", ...
            "The field %s must not be empty.", field_name);
    end
end
