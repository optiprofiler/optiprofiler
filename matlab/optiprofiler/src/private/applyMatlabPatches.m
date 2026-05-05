function guard = applyMatlabPatches()
%APPLYMATLABPATCHES  Temporary, self-restoring workarounds for known
%   MATLAB bugs that affect OptiProfiler at runtime.
%
%   This is a single, narrow place to collect small "patches" for issues
%   that live in MATLAB itself (or in the surrounding system) rather than
%   in OptiProfiler. Each patch:
%     * detects its own risky environment and otherwise returns early,
%     * applies the smallest possible runtime change,
%     * registers a restore callback so the original state is recovered
%       on exit.
%
%   The function returns a single onCleanup handle. When that handle goes
%   out of scope (i.e. when the calling function returns, normally or via
%   an error), every registered restore callback is invoked in reverse
%   order. No MATLAB preferences are written and no persistent state is
%   modified, so this helper is safe to call unconditionally.
%
%   Usage at the top of any entry point:
%
%       patchesGuard = applyMatlabPatches();  %#ok<NASGU>
%
%   The guard MUST remain a local variable so it is destroyed on return.
%
%   Adding a new patch:
%     1. Implement a `local<Name>` sub-function that returns its updated
%        cell array of restorers (or the input cell unchanged if the
%        environment is not affected).
%     2. Add one line to the main function below:
%            restorers = local<Name>(restorers);
%
%   To remove a patch, simply delete its sub-function and the matching
%   line. The helper is intentionally easy to grow and easy to retire.

    restorers = {};

    % Patch: MATLAB R2025b on Linux with software OpenGL fallback can
    % crash inside the new web-based graphics backend (libmwhgweb,
    % `hg::web::BinaryIdToCallbackRouter` on the SceneTree thread).
    restorers = localPatchR2025bLinuxGraphics(restorers);

    % --- Add future patches here, e.g.: -----------------------------------
    % restorers = localPatchSomeOtherIssue(restorers);

    guard = onCleanup(@() localRunRestorers(restorers));
end


function restorers = localPatchR2025bLinuxGraphics(restorers)
% Mitigation for an intermittent crash in MATLAB R2025b on Linux when
% MATLAB has fallen back to software OpenGL. The crash is an internal
% MATLAB assertion in the new web-based graphics backend (libmwhgweb,
% `hg::web::BinaryIdToCallbackRouter` on the SceneTree thread) and is
% triggered by the desktop figure infrastructure itself, not by any
% specific OptiProfiler call. Switching the default renderer to
% 'painters' and rendering off-screen reduces the likelihood; it does
% not eliminate it. The only guaranteed workaround is to run MATLAB
% without the desktop GUI (e.g. via `matlab -batch ...`).

    if ~isunix() || ismac()
        return;
    end
    if exist('isMATLABReleaseOlderThan', 'file') ~= 2
        return;
    end
    if isMATLABReleaseOlderThan('R2025b')
        return;
    end
    % Scope this workaround to R2025b only: newer yearly releases may fix
    % hgweb; if reports show the bug persists, widen this bound explicitly.
    if ~isMATLABReleaseOlderThan('R2026a')
        return;
    end
    try
        d = opengl('data');
        software_opengl = isstruct(d) && isfield(d, 'Software') ...
            && logical(d.Software);
    catch
        software_opengl = false;
    end
    if ~software_opengl
        return;
    end

    g = groot;
    saved.Renderer     = get(g, 'defaultFigureRenderer');
    saved.RendererMode = get(g, 'defaultFigureRendererMode');
    saved.Visible      = get(g, 'defaultFigureVisible');

    set(g, ...
        'defaultFigureRenderer',     'painters', ...
        'defaultFigureRendererMode', 'manual', ...
        'defaultFigureVisible',      'off');

    warning('OptiProfiler:MatlabPatchApplied', ...
        ['Detected MATLAB R2025b on Linux with software OpenGL. ' ...
         'In this configuration MATLAB itself can crash intermittently ' ...
         'during desktop figure rendering (an internal assertion in ' ...
         'libmwhgweb / hg::web::BinaryIdToCallbackRouter on the SceneTree ' ...
         'thread). This is a MATLAB issue, independent of OptiProfiler. ' ...
         'As a best-effort mitigation, the default figure renderer has ' ...
         'been temporarily set to ''painters'' and new figures forced to ' ...
         'invisible; the original defaults will be restored automatically ' ...
         'on exit. The mitigation reduces the likelihood but cannot ' ...
         'eliminate it. The only guaranteed workaround is to run MATLAB ' ...
         'without the desktop GUI, i.e. in headless / batch mode.']);

    restorers{end + 1} = @() localRestoreGrootDefaults(saved);
end


function localRestoreGrootDefaults(saved)
    try
        g = groot;
        set(g, ...
            'defaultFigureRenderer',     saved.Renderer, ...
            'defaultFigureRendererMode', saved.RendererMode, ...
            'defaultFigureVisible',      saved.Visible);
    catch
        % Never propagate errors out of cleanup.
    end
end


function localRunRestorers(restorers)
    for i = numel(restorers):-1:1
        try
            restorers{i}();
        catch
            % Never propagate errors out of cleanup.
        end
    end
end
