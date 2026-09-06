function [unsafe, reason] = nativeGraphicsTempPathIsUnsafe(tmpdir_value, platform, release, tmpdir_is_set)
%NATIVEGRAPHICSTEMPPATHISUNSAFE Recognize one verified native-renderer risk.
%   R2026a on Linux appends the 46-byte ASCII path component
%   '/.org.chromium.Chromium.XXXXXX/SingletonSocket' to TMPDIR. Linux's
%   sockaddr_un.sun_path[108] leaves at most 107 bytes before the final NUL.
%   The affected child can terminate while exportgraphics keeps waiting.
%   This is not a universal tempdir limit: preserve other releases/platforms.
%   Inputs are explicit so this policy can be tested without changing globals.
    unsafe = false;
    reason = '';
    if ~strcmp(platform, 'glnxa64') || ~strcmp(release, '2026a')
        return;
    end
    if isempty(tmpdir_value)
        % MATLAB getenv returns '' both for an unset variable (safe default)
        % and an explicitly empty value (a reproduced renderer failure).
        % The caller supplies live presence, or [] if that check failed.
        if nargin < 4 || isempty(tmpdir_is_set)
            unsafe = true;
            reason = 'TMPDIR is empty and its environment presence could not be established.';
        elseif tmpdir_is_set
            unsafe = true;
            reason = 'TMPDIR is explicitly set to an empty value, which the native renderer cannot safely use.';
        end
        return;
    end
    tmpdir_value = char(tmpdir_value);
    if tmpdir_value(1) ~= '/'
        % Conservative policy, not a claim that every relative path hangs:
        % the native child need not resolve a relative path like MATLAB does.
        unsafe = true;
        reason = 'TMPDIR is relative, so the native renderer path cannot be safely established.';
        return;
    end
    % Appending a component removes redundant trailing separators. Keep '/'.
    tmpdir_value = regexprep(tmpdir_value, '/+$', '');
    if isempty(tmpdir_value), tmpdir_value = '/'; end
    byte_count = numel(unicode2native(tmpdir_value, 'UTF-8'));
    socket_suffix_bytes = 46;
    socket_path_capacity = 107;
    unsafe = byte_count + socket_suffix_bytes > socket_path_capacity;
    if unsafe
        reason = sprintf('TMPDIR uses %d UTF-8 bytes; the observed renderer suffix exceeds the Linux socket path capacity.', byte_count);
    end
end
