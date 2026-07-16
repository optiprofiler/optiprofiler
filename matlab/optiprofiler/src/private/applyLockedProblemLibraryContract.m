function registration = applyLockedProblemLibraryContract(registration)
%APPLYLOCKEDPROBLEMLIBRARYCONTRACT applies canonical metadata to known names.

    [lock, ~] = readMatlabProblemLibraryLock();
    entries = lock.libraries;
    for i_entry = 1:numel(entries)
        if iscell(entries)
            entry = entries{i_entry};
        else
            entry = entries(i_entry);
        end
        if ~strcmp(char(entry.name), registration.name)
            continue
        end
        if strcmp(char(entry.role), 'bundled-default')
            error("MATLAB:problemLibraryRegistry:bundledRegistration", ...
                "The bundled default '%s' is resolved from its locked gitlink.", ...
                registration.name);
        end
        if ~strcmp(char(entry.select_function), registration.select_function) || ...
                ~strcmp(char(entry.load_function), registration.load_function)
            error("MATLAB:problemLibraryRegistry:lockedContractMismatch", ...
                "Registration for '%s' does not use its locked select/load functions.", ...
                registration.name);
        end
        registration.api_version = double(entry.api_version);
        registration.role = char(entry.role);
        registration.collect_info_function = optionalString(entry.collect_info_function);
        registration.platforms = cellstr(entry.platforms);
        registration = normalizeProblemLibraryRegistration(registration);
        return
    end
end


function value = optionalString(value)
    if isempty(value)
        value = '';
    else
        value = char(value);
    end
end
