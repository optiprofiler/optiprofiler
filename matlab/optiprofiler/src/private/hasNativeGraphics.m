function available = hasNativeGraphics()
%HASNATIVEGRAPHICS Probe graphics capability independently of JVM presence.
%   Cached per MATLAB process: do not repeatedly allocate a probe figure for
%   every history plot. No Java, persistent path change or external program.
    persistent result
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
