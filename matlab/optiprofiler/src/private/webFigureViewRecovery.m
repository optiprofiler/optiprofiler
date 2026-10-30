function status = webFigureViewRecovery(fig)
%WEBFIGUREVIEWRECOVERY reloads the page of a figure whose web view is not ready.
%   STATUS = WEBFIGUREVIEWRECOVERY(FIG) is called once by boundedFigureCall while
%   a native export of FIG has been waiting for some time. It reloads FIG's page
%   if all of the following hold, and otherwise changes nothing:
%   - MATLAB R2026a on Linux (the release and platform where recovery passed
%     validation; R2025a recovery is withheld after unresolved heap errors);
%   - FIG is a valid invisible figure (a visible one may be a window the user
%     works in) whose view is still not ready;
%   - exactly one window of MATLAB's embedded browser shows FIG's own page from
%     MATLAB's local web server.
%   That window is pointed at the same page again, without the expired one-time
%   access parameter, to which MATLAB's window object adds a fresh one. The
%   function then waits until the reloaded view is ready (at most 20 seconds),
%   and 1 second more, so that the export resumes with a settled layout. Nothing
%   is closed or deleted, and no other window is touched. STATUS names the
%   outcome. Errors are contained, and an unavailable or unmatched window is
%   left untouched. The caller's cooperative time limit remains in force.
%
%   A new figure's view becomes ready only after its page has loaded from
%   MATLAB's local web server. If one of the page's own requests fails, the page
%   never reports ready, and an export waiting for the view waits forever (see
%   boundedFigureCall). Loading the page again normally completes the view, and
%   the waiting export continues. A view that is merely slow is also reloaded,
%   which costs one extra page load. Because this relies on undocumented MATLAB
%   internals (the figure's Uuid and FigureViewReady properties, the window list
%   of matlab.internal.cef.webwindowmanager and the address of the figure page),
%   it is limited to the releases and platform where it was tested. Recheck it
%   with the TestRuntimeRobustness tests before adding a release.
    status = 'unsupported';
    try
        if ~webFigureRecoverySupported()
            return;
        end
        if ~isscalar(fig) || ~isgraphics(fig, 'figure')
            status = 'figure gone';
            return;
        end
        if ~strcmp(char(string(fig.Visible)), 'off')
            return;
        end
        if isReady(fig)
            status = 'view ready';
            return;
        end
        % The figure page is .../cefComponentContainer.html?channel=/uifigure/<Uuid>&...
        % on MATLAB's local web server; nothing else is ever reloaded.
        page = ['^https?://(127\.0\.0\.1|localhost):\d+/toolbox/matlab/uitools/uifigureappjs/' ...
            'cefComponentContainer\.html\?(.*&)?channel=/uifigure/', ...
            regexptranslate('escape', char(string(fig.Uuid))), '(&|$)'];
        windows = matlab.internal.cef.webwindowmanager.instance.windowList;
        matches = false(size(windows));
        for k = 1:numel(windows)
            matches(k) = ~isempty(regexp(char(windows(k).URL), page, 'once'));
        end
        if nnz(matches) ~= 1
            status = 'no unique figure page';
            return;
        end
        window = windows(matches);
        if ~window.isWindowValid
            status = 'window closed';
            return;
        end
        url = regexprep(char(window.URL), '([?&])snc=[^&#]*&?', '$1');
        window.URL = regexprep(url, '[?&]$', '');
        status = 'reloaded';
        settle(fig);
    catch
        status = 'not reloaded';
    end
end

function settle(fig)
% Waits until the reloaded view is ready, and then briefly, before the waiting
% export resumes. While the page was broken, MATLAB laid out the figure with
% provisional values that the page corrects only once it runs; an export
% resuming at once captured them (a log-ratio profile came out about 5% shorter).
% If the view is still not ready after 20 seconds, it returns at once.
    try
        waited = tic;
        while toc(waited) < 20 && isgraphics(fig) && ~isReady(fig)
            pause(0.05);
        end
        if isgraphics(fig) && isReady(fig)
            drawnow;
            pause(1);
        end
    catch
    end
end

function ready = isReady(fig)
    ready = strcmp(char(string(get(fig, 'FigureViewReady'))), 'on');
end
