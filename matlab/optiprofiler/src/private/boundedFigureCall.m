function boundedFigureCall(fig, call, limit)
%BOUNDEDFIGURECALL runs CALL, a native export of FIG, with a time limit.
%   In MATLAB's web-based graphics, a native export of a new figure first waits,
%   without a timeout, until the figure's page in the embedded Chromium browser
%   reports that its view is ready. If one of the page's own requests to MATLAB's
%   local web server fails, the page never reports ready and the export waits
%   indefinitely. This was observed on Linux with R2025a and R2026a, also in plain
%   MATLAB code, when Chromium cancelled such a request because the host's network
%   configuration changed. MATLAB keeps running timer callbacks during that wait,
%   so a watchdog timer acts on it:
%   - On Linux R2026a only, after RECOVERY_DELAY seconds, if CALL is still running, it calls
%     webFigureViewRecovery once, which reloads FIG's page where that is
%     supported and FIG's view is still not ready, and waits for the reloaded
%     view; the waiting export then continues with it.
%   - After LIMIT seconds (default 300, about a hundred times an ordinary
%     export), it deletes FIG. That ends the wait, and the call fails with
%     OptiProfiler:NativeExportTimeout.
%   A healthy export normally finishes long before the first event, and a slower
%   one whose view is ready is not reloaded. MATLAB timers are cooperative: they
%   act only while MATLAB processes events, as it did during the observed waits
%   for page initialization, and cannot bound a call that blocks MATLAB's event loop.
%   Without the watchdog, an export could wait forever, so CALL is not run when
%   the watchdog cannot be started (OptiProfiler:NativeExportUnbounded). Callers
%   turn native-export errors into the portable SVG output. Errors of CALL
%   itself are rethrown unchanged.
%   Tests: TestRuntimeRobustness (bounded native export tests).
    if nargin < 3
        limit = 300;
    end
    % Healthy figure views were ready well within this delay in all
    % measurements; a view that is still not ready then is reloaded once.
    RECOVERY_DELAY = 30;
    watchdog = [];
    try
        % The callback captures only the figure: a nested function would keep
        % this workspace, including the cleanup below, alive with the timer.
        if webFigureRecoverySupported() && limit > RECOVERY_DELAY
            % Two events: the reload attempt, then the time limit.
            watchdog = timer('Tag', 'OptiProfilerNativeExportWatchdog', 'ExecutionMode', 'fixedRate', ...
                'StartDelay', RECOVERY_DELAY, 'Period', limit - RECOVERY_DELAY, 'TasksToExecute', 2, ...
                'UserData', struct('reload', true, 'expired', false), 'TimerFcn', @(source, ~) act(source, fig));
        else
            watchdog = timer('Tag', 'OptiProfilerNativeExportWatchdog', 'ExecutionMode', 'singleShot', ...
                'StartDelay', limit, 'UserData', struct('reload', false, 'expired', false), ...
                'TimerFcn', @(source, ~) act(source, fig));
        end
        start(watchdog);
    catch cause
        removeWatchdog(watchdog);
        error('OptiProfiler:NativeExportUnbounded', ['The native export was not attempted because no ' ...
            'MATLAB timer could limit its duration: %s'], cause.message);
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
        error('OptiProfiler:NativeExportTimeout', ...
            'The native export did not finish within %g seconds, so the figure was closed.', limit);
    end
end

function act(watchdog, fig)
    state = watchdog.UserData;
    if state.reload
        state.reload = false;
        watchdog.UserData = state;
        % An error in a timer callback stops the timer, which would also
        % cancel the time limit, so no failure of the reload may escape. The
        % reload waits at most 20 seconds for the view, long before the limit.
        try
            webFigureViewRecovery(fig);
        catch
        end
        return;
    end
    state.expired = true;
    watchdog.UserData = state;
    if isgraphics(fig)
        delete(fig);
    end
end

function expired = hasExpired(watchdog)
    expired = ~isempty(watchdog) && isvalid(watchdog) && watchdog.UserData.expired;
end

function removeWatchdog(watchdog)
    if ~isempty(watchdog) && isvalid(watchdog)
        stop(watchdog);
        delete(watchdog);
    end
end
