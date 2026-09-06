function atomicReplaceFile(source, target)
%ATOMICREPLACEFILE Publish a completed file without deleting the previous one.
%   Java uses ATOMIC_MOVE. On macOS/Linux without Java, mktemp owns a staging
%   file in the target's SAME directory; mv can then use the POSIX rename
%   primitive, never its cross-filesystem copy-and-delete fallback.
    if ~isfile(source) || isfolder(target) || strcmp(source,target)
        error('OptiProfiler:AtomicReplace', 'Source must be a file and target must not be a directory.');
    end
    if usejava('jvm')
        from = java.io.File(source); to = java.io.File(target);
        options = javaArray('java.nio.file.CopyOption', 2);
        options(1) = java.nio.file.StandardCopyOption.ATOMIC_MOVE;
        options(2) = java.nio.file.StandardCopyOption.REPLACE_EXISTING;
        java.nio.file.Files.move(from.toPath(), to.toPath(), options);
    elseif isunix
        if ~startsWith(source, filesep), source = fullfile(pwd, source); end
        if ~startsWith(target, filesep), target = fullfile(pwd, target); end
        [source_ok, source_parent] = fileattrib(fileparts(source));
        [target_ok, target_parent] = fileattrib(fileparts(target));
        same_parent = source_ok && target_ok && strcmp(source_parent.Name,target_parent.Name);
        staging = source;
        if ~same_parent
            % Raw MAT staging already shares the target directory: avoid a
            % redundant multi-GB copy. Only cross-directory publication needs
            % this extra same-parent staging file to exclude mv's copy mode.
            template = fullfile(fileparts(target), '.optiprofiler-replace.XXXXXXXX');
            [status, output] = system(['/usr/bin/mktemp ', quotePosix(template), ' 2>&1']);
            if status ~= 0
                error('OptiProfiler:AtomicReplace', 'Cannot reserve replacement file: %s', strtrim(output));
            end
            staging = strtrim(output);
            try
                copyfile(source, staging, 'f');
            catch cause
                removeStaging(staging); % Only incomplete copies are disposable.
                rethrow(cause);
            end
        end
        [status, output] = system(['/bin/mv -f ', quotePosix(staging), ' ', quotePosix(target), ' 2>&1']);
        if status ~= 0
            % A complete copy is recoverable even when publishing fails.
            error('OptiProfiler:AtomicReplace', 'Cannot replace output file: %s. Completed staging retained at ''%s''.', strtrim(output), staging);
        end
        if ~same_parent, delete(source); end
    else
        error('OptiProfiler:SafeOutputPlatform', 'Enable the JVM for atomic output replacement on this platform.');
    end
end

function removeStaging(file)
    if isfile(file), delete(file); end
end
