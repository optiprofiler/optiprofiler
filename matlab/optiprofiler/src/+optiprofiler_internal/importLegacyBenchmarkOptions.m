function [feature, receipt] = importLegacyBenchmarkOptions(source, legacy_feature_name)
%IMPORTLEGACYBENCHMARKOPTIONS Convert known unversioned native settings.
% The public reader owns files and replay options; this adapter owns old feature
% layouts and retains only facts that the saved configuration actually records.
    receipt = struct();
    if isfield(source, 'feature')
        if isa(source.feature, 'optiprofiler_internal.LegacyFeatureEnvelope')
            [feature, request] = optiprofiler_internal.importLegacyFeature(source.feature);
            receipt.feature_import = request;
        else
            feature = Feature(source.feature);
        end
    else
        if isfield(source, 'feature_name')
            name = source.feature_name;
        elseif ~isempty(legacy_feature_name)
            name = legacy_feature_name;
        else
            error('OptiProfiler:LegacyFeatureIdentityMissing', ...
                'This old flat options file did not retain feature identity. Supply legacy_feature_name explicitly; stamps and folders are not evidence.');
        end
        definitions = optiprofiler_internal.featureDefinitions();
        local_keys = {};
        for i = 1:numel(definitions), local_keys = [local_keys, definitions(i).local_keys]; end %#ok<AGROW>
        local_options = struct();
        keys = intersect(fieldnames(source), unique(local_keys));
        for i = 1:numel(keys), local_options.(keys{i}) = source.(keys{i}); end
        saved = struct('name', name, 'options', local_options);
        if isfield(source, 'n_runs'), saved.options.n_runs = source.n_runs; end
        [feature, request] = optiprofiler_internal.importLegacyFeature(optiprofiler_internal.LegacyFeatureEnvelope(saved));
        receipt.feature_import = request;
        if ~isempty(legacy_feature_name) && ~isfield(source, 'feature_name')
            receipt.current_legacy_feature_name_override = legacy_feature_name;
        end
    end
end
