function boundedFigureCall(fig, call, limit)
%BOUNDEDFIGURECALL runs CALL, a native export of FIG, for at most LIMIT seconds.
%   MATLAB's web-based graphics can wait indefinitely for its renderer during
%   an export (observed on Linux with R2025a and R2026a, also in plain MATLAB
%   code without OptiProfiler). MATLAB keeps running timer callbacks during
%   that wait, and deleting the figure ends it. A one-shot timer therefore
%   deletes FIG when CALL has not returned after LIMIT seconds (default 300,
%   about a hundred times an ordinary export), and the call then fails with
%   OptiProfiler:NativeExportTimeout. Callers already turn native-export
%   errors into the portable SVG output, so results are never lost.
%   Errors of CALL itself are rethrown unchanged.
    if nargin < 3
        limit = 300;
    end
    watchdog = [];
    try
        % The callback captures only the figure: a nested function would keep
        % this workspace, including the cleanup below, alive with the timer.
        watchdog = timer('Tag', 'OptiProfilerNativeExportWatchdog', 'ExecutionMode', 'singleShot', ...
            'StartDelay', limit, 'UserData', false, 'TimerFcn', @(source, ~) expire(source, fig));
        start(watchdog);
    catch
        % Timers are unavailable (for example, in a thread-based worker):
        % export without a bound, as before.
        removeWatchdog(watchdog);
        watchdog = [];
    end
    cleanup = onCleanup(@() removeWatchdog(watchdog));
    try
        call();
    catch cause
        if ~hasExpired(watchdog)
            rethrow(cause);
        end
    end
    if hasExpired(watchdog)
        error('OptiProfiler:NativeExportTimeout', ['The native export did not finish within %g seconds; ' ...
            'MATLAB''s graphics renderer stopped responding, so the figure was closed.'], limit);
    end
end

function expire(watchdog, fig)
    watchdog.UserData = true;
    if isgraphics(fig)
        delete(fig);
    end
end

function expired = hasExpired(watchdog)
    expired = ~isempty(watchdog) && isvalid(watchdog) && isequal(watchdog.UserData, true);
end

function removeWatchdog(watchdog)
    if ~isempty(watchdog) && isvalid(watchdog)
        stop(watchdog);
        delete(watchdog);
    end
end
