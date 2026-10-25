function appendReferenceHistoryDiagnostics(path_log, result, library, role)
%APPENDREFERENCEHISTORYDIAGNOSTICS Persist silent-trial causes on the controller.
% Workers only return receipts. This runs once after collection, never writes
% to the console, and does not duplicate copied runs or ordinary-mode warnings.
    lines = diagnosticLines(result, library, role);
    if isempty(lines), return; end
    target = fullfile(path_log, 'log.txt');
    [stream, message] = fopen(target, 'a');
    if stream < 0
        error('OptiProfiler:ReferenceHistoryDiagnosticWrite', ...
            'Cannot append reference-history diagnostics: %s', message);
    end
    cleanup = onCleanup(@() fclose(stream)); %#ok<NASGU>
    for k = 1:numel(lines)
        fprintf(stream, '%s\n', lines{k});
    end
end

function lines = diagnosticLines(result, library, role)
    lines = {};
    if isfield(result, 'eval_report_metadata')
        lines = trialLines(result.eval_report_metadata, result.problem_name, library, role);
    elseif isfield(result, 'execution_metadata')
        for p = 1:numel(result.execution_metadata)
            lines = [lines, trialLines(result.execution_metadata{p}, ...
                result.problem_names{p}, library, role)]; %#ok<AGROW>
        end
    end
    if isfield(result, 'results_plib_plain')
        lines = [lines, diagnosticLines(result.results_plib_plain, library, 'plain_reference')];
    end
end

function lines = trialLines(metadata, problem_name, library, role)
    lines = {};
    if ~isfield(metadata, 'runtime_receipts'), return; end
    receipts = metadata.runtime_receipts;
    for solver = 1:size(receipts, 1)
        for run = 1:size(receipts, 2)
            receipt = receipts{solver, run};
            if ~isstruct(receipt) || ~isfield(receipt, 'reference_history_failure'), continue; end
            failure = receipt.reference_history_failure;
            lines{end+1} = sprintf( ...
                '[WARNING] Reference history unavailable for %s/%s (%s, solver %d, run %d): %s: %s', ...
                library, shortenMessageForLog(problem_name), role, solver, run, ...
                failure.identifier, failure.message); %#ok<AGROW>
        end
    end
end
