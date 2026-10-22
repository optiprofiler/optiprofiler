function supported = webFigureRecoverySupported()
%WEBFIGURERECOVERYSUPPORTED scopes the internal figure-page repair.
%   Keep both the watchdog schedule and the repair itself behind this gate.
%   R2025a reloads had unresolved post-export heap errors; other releases have
%   not been validated. They retain the ordinary single-shot timeout instead.
    supported = strcmp(computer('arch'), 'glnxa64') && strcmp(version('-release'), '2026a') ...
        && exist('matlab.internal.cef.webwindowmanager', 'class') == 8;
end
