function [feature, experiment_request] = importLegacyFeature(envelope)
%IMPORTLEGACYFEATURE Explicit trusted import of a frozen ac4 native Feature.
% [FEATURE, REQUEST] = optiprofiler_internal.importLegacyFeature(ENVELOPE)
% separates the retained resolved n_runs from clean local stages. REQUEST does
% not assert that this value was explicitly requested by the historical user.
% Unrecorded declarations remain unknown. JSON callback labels are not inputs.
    if ~strcmp(class(envelope),'optiprofiler_internal.LegacyFeatureEnvelope') || ~isscalar(envelope)
        error('MATLAB:Feature:LegacyImportRequired','Import requires the legacy envelope returned by trusted native load.');
    end
    saved = envelope.saved_state;
    if ~(isa(saved.name,'char') && isrow(saved.name)) && ~(isa(saved.name,'string') && isscalar(saved.name))
        error('MATLAB:Feature:InvalidNativeState','Historical feature identity must be present as scalar text.');
    end
    name = lower(strtrim(char(saved.name)));
    if isempty(name) || contains(name,'+') || ~isa(saved.options,'struct') || ~isscalar(saved.options)
        error('MATLAB:Feature:InvalidNativeState','Only the known atomic name/options historical layout is supported.');
    end
    local = saved.options;
    experiment_request = struct('schema','matlab-legacy-experiment-request-v1', ...
        'n_runs',[],'n_runs_present',false,'origin','unknown','requested_n_runs',[], ...
        'requested_origin','unknown','source_format','ac4-feature-name-options');
    if isfield(local,'n_runs')
        value = local.n_runs;
        if ~builtin('isnumeric',value) || ~isreal(value) || ~isscalar(value) || rem(value,1)~=0 || value<=0
            error('MATLAB:Feature:InvalidNativeState','The retained resolved n_runs is not a positive integer.');
        end
        experiment_request.n_runs = value;
        experiment_request.n_runs_present = true;
        experiment_request.origin = 'retained_resolved_value';
        local = rmfield(local,'n_runs');
    end
    definition = optiprofiler_internal.featureDefinitions(name);
    if strcmp(name,'plain')
        if ~isempty(fieldnames(local))
            error('MATLAB:Feature:InvalidNativeState','Historical plain Feature has unsupported local options.');
        end
        stages = {};
    else
        stages = {struct('name',name,'options',local,'occurrence',0, ...
            'identity',[name,'#0'],'code',definition.code)};
    end
    state = struct('stages',{stages},'declared',[],'declared_name','', ...
        'specification_version','matlab-feature-spec-v2');
    % The current native boundary performs the one validation/normalization
    % pass; the import adapter does not construct-and-discard a Feature.
    feature = Feature.loadobj(struct('schema','matlab-feature-native-v2','specification',state));
end
