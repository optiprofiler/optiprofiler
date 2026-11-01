function directory = resolvedDirectory(directory)
%RESOLVEDDIRECTORY returns the name that the MATLAB path holds for a directory.
%   The directory may not exist yet: what does exist of it is resolved.

    missing = {};
    existing = directory;
    while ~isfolder(existing)
        [parent, name, extension] = fileparts(existing);
        if strcmp(parent, existing) || isempty(parent)
            return
        end
        missing = [{[name, extension]}, missing]; %#ok<AGROW>
        existing = parent;
    end
    if isunix
        % addpath follows a link in the last component, fileattrib does not.
        existing = fullfile(existing, '.');
    end
    [status, attributes] = fileattrib(existing);
    if status
        directory = fullfile(attributes.Name, missing{:});
    end
end
