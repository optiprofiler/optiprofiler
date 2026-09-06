function available = hasNativeGraphics()
%HASNATIVEGRAPHICS Probe graphics capability independently of JVM presence.
%   Cached per MATLAB process: do not repeatedly allocate a probe figure for
%   every history plot. No Java or environment/path changes. Only an empty
%   Linux R2026a TMPDIR needs a child presence check (getenv is ambiguous).
    persistent result warned_unsafe_tmpdir
    % Recheck before the figure cache: a caller may change TMPDIR later.
    % A policy rejection must not poison the real figure-capability cache.
    tmpdir_value = getenv('TMPDIR');
    platform = computer('arch');
    release = version('-release');
    tmpdir_is_set = [];
    if strcmp(platform, 'glnxa64') && strcmp(release, '2026a') && isempty(tmpdir_value)
        % Constant child command: inspect only presence, never print the
        % environment or change it. Query each time; Java/proc snapshots can
        % be stale after user callbacks call setenv. Status 0=set, 1=unset.
        try
            [status, ~] = system('/bin/sh -c ''test "${TMPDIR+x}" = x''');
            if status == 0, tmpdir_is_set = true;
            elseif status == 1, tmpdir_is_set = false;
            end
        catch
            % An unknown presence stays [] and takes conservative SVG output.
        end
    end
    [unsafe, reason] = nativeGraphicsTempPathIsUnsafe(tmpdir_value, platform, release, tmpdir_is_set);
    if unsafe
        available = false;
        if isempty(warned_unsafe_tmpdir)
            warned_unsafe_tmpdir = true;
            % Routing must remain visible even while history plotting has
            % warnings disabled, and must not become an error via warning policy.
            printOptiProfilerMessage('WARNING', sprintf( ...
                ['%s Avoiding a known MATLAB R2026a Linux renderer hang; portable SVG/HTML output remains available. ', ...
                'Use a short absolute TMPDIR before plotting or starting MATLAB. OptiProfiler has not changed TMPDIR.'], reason));
        end
        return;
    end
    if isempty(result)
        try
            f = figure('Visible', 'off');
            close(f);
            result = true;
        catch
            result = false;
        end
    end
    available = result;
end
