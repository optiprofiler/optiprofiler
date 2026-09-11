function records = retainedExecutionMetadata(results)
%RETAINEDEXECUTIONMETADATA Keep observations, never infer calls from a plan.
% One cell per retained problem. Do not serialize transient graphics, handles,
% raw callbacks or the report's presentation caches into numerical archives.
    records = cell(size(results));
    allowed = {'real_n_runs', 'oracle_seeds', 'runtime_receipts'};
    for p = 1:numel(results)
        record = struct();
        if isfield(results{p}, 'eval_report_metadata')
            metadata = results{p}.eval_report_metadata;
            for k = 1:numel(allowed)
                key = allowed{k};
                if isfield(metadata, key), record.(key) = metadata.(key); end
            end
        end
        records{p} = record;
    end
end
