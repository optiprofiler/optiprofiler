function platform = matlabPlatformName()
%MATLABPLATFORMNAME returns the registry platform identifier.

    if ispc
        platform = 'windows';
    elseif ismac
        platform = 'macos';
    else
        platform = 'linux';
    end
end
