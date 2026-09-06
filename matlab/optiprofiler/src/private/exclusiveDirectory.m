function created = exclusiveDirectory(directory)
%EXCLUSIVEDIRECTORY Acquire a directory exactly once; never claim an old one.
    if usejava('jvm')
        object = java.io.File(directory);
        created = object.mkdir();
        exists = object.exists();
        message = 'Java File.mkdir failed.';
    elseif isunix
        % POSIX mkdir is the ownership operation; MATLAB mkdir reports success
        % for an existing directory and therefore cannot reserve output safely.
        if ~startsWith(directory, filesep), directory = fullfile(pwd, directory); end
        [status, message] = system(['/bin/mkdir ', quotePosix(directory), ' 2>&1']);
        created = status == 0;
        exists = isfolder(directory) || isfile(directory);
    else
        error('OptiProfiler:SafeOutputPlatform', ...
            'No-JVM safe output is supported on macOS/Linux. Enable the JVM on this platform.');
    end
    if ~created && ~exists
        error('OptiProfiler:OutputDirectory', 'Cannot create directory ''%s'': %s', directory, strtrim(message));
    end
end
