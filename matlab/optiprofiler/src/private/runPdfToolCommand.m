function [status, output, retried] = runPdfToolCommand(command, temp_output)
%RUNPDFTOOLCOMMAND Retry only an identified Linux C++ loader incompatibility.
%   TEMP_OUTPUT is the caller's private staging file, never the published
%   summary or an input. Tool order, arguments and publication stay with the
%   caller. No MATLAB environment variable is changed, even while retrying.

    persistent announced_retry;
    command = string(command);
    [status, output] = system(command + " 2>&1");
    output = string(strtrim(output));
    retried = false;
    % MATLAB's bundled libstdc++ may be older than the OS PDF tool requires.
    % Match the loader's specific version diagnostic, not PDF syntax errors,
    % a missing executable, or arbitrary mentions of a version in document text.
    loader_failure = "(?m)^(?:[^:\r\n]+:\s*)?/[^\r\n]*libstdc\+\+\.so(?:\.[0-9]+)*:\s*" + ...
        "version [`'](?:GLIBCXX|CXXABI)_[0-9]+(?:\.[0-9]+)*['`]\s+not found \(required by [^\r\n]+\)";
    if status == 0 || ~isunix || ismac || isempty(getenv('LD_LIBRARY_PATH')) || ...
            isempty(regexp(char(output), char(loader_failure), 'once'))
        return;
    end

    % A wrapper can create a partial file before exec reaches the loader.
    % Remove only this attempt's staging file: exit-zero with no new output
    % must not cause the caller to publish that first, failed partial copy.
    if isfile(temp_output)
        try
            delete(temp_output);
        catch cause
            output = output + newline + "Cannot remove the failed temporary PDF before loader retry: " + string(cause.message);
            return;
        end
        if isfile(temp_output)
            output = output + newline + "Cannot remove the failed temporary PDF before loader retry: " + string(temp_output);
            return;
        end
    end

    % Quote the whole existing command once for the clean child shell. In
    % particular, quoted PDF filenames must not become shell substitutions.
    quoted = "'" + replace(command + " 2>&1", "'", "'\''") + "'";
    [status, retry_output] = system("/usr/bin/env -u LD_LIBRARY_PATH /bin/sh -c " + quoted + " 2>&1");
    retried = true;
    if status == 0
        output = string(strtrim(retry_output));
        if isfile(temp_output) && isempty(announced_retry)
            printOptiProfilerMessage('INFO', ...
                'PDF tool subprocess retry without LD_LIBRARY_PATH succeeded; parent MATLAB environment unchanged.');
            announced_retry = true;
        end
    else
        output = "Inherited environment: " + output + newline + ...
            "Retry without LD_LIBRARY_PATH failed: " + string(strtrim(retry_output));
    end
end
