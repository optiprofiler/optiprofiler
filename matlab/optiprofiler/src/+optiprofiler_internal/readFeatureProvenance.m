function [payload, interpretation] = readFeatureProvenance(value)
%READFEATUREPROVENANCE Retain historical payloads; never upgrade on read.
% Unsupported/malformed metadata does not require re-running a provider or
% solver. Numeric load can proceed while the interpretation remains explicit.
    unknown = struct('eval_report_null', true);
    interpretation = struct('status', 'absent', 'schema', unknown);
    payload = unknown;
    if nargin == 0, return; end
    payload = value;
    if ischar(value) || (isstring(value) && isscalar(value))
        try
            payload = jsondecode(char(value));
        catch
            interpretation.status = 'malformed';
            interpretation.reason = 'unparsable_feature_pipeline';
            return;
        end
    end
    if ~isstruct(payload) || ~isscalar(payload) || ~isfield(payload, 'schema') ...
            || ~(ischar(payload.schema) || (isstring(payload.schema) && isscalar(payload.schema)))
        interpretation.status = 'malformed';
        interpretation.reason = 'missing_feature_pipeline_schema';
        return;
    end
    interpretation.schema = payload.schema;
    if ismember(payload.schema, {'feature_pipeline-v1', 'feature_pipeline-v2', 'feature_pipeline-v3'})
        interpretation.status = 'known';
    else
        interpretation.status = 'unsupported';
        interpretation.reason = 'unknown_feature_pipeline_schema';
    end
end
