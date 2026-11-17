function [name, vouched] = archivedFeatureName(result)
%ARCHIVEDFEATURENAME returns the name of the feature of an archived library record, for display.
%   It is the effective name recorded in the `feature_pipeline` entry (any known version). An
%   archive without that entry (OptiProfiler 1.x) only has its `feature_stamp`: a generated stamp is
%   the feature name followed by option values ('linearly_transformed_rotated'; the stages of a
%   composition are joined with a double underscore), so the name is read from it. A stamp of
%   another form (it can be set by the user) is returned as it is, since it is the only name the
%   archive has. The name is empty when the record has neither.
%
%   `vouched` lists the feature names that the record vouches for, to compare a given name with:
%   the effective and declared names of the `feature_pipeline` entry, otherwise the name read from
%   a generated stamp. It is empty when the record has neither, or only a stamp of another form.

    name = '';
    vouched = {};
    if isfield(result, 'feature_pipeline')
        payload = optiprofiler_internal.readFeatureProvenance(result.feature_pipeline);
        if isstruct(payload) && isscalar(payload)
            % The feature block of the current schema; the first schema was flat.
            block = payload;
            if isfield(payload, 'feature') && isstruct(payload.feature) && isscalar(payload.feature)
                block = payload.feature;
            end
            if isfield(block, 'effective_name') && isScalarText(block.effective_name) && strlength(block.effective_name) > 0
                name = char(block.effective_name);
                vouched = {name};
                if isfield(block, 'declared_name') && isScalarText(block.declared_name) && strlength(block.declared_name) > 0
                    vouched = unique([vouched, {char(block.declared_name)}], 'stable');
                end
                return;
            end
        end
    end
    if ~isfield(result, 'feature_stamp') || ~isScalarText(result.feature_stamp) || strlength(result.feature_stamp) == 0
        return;
    end
    name = char(result.feature_stamp);
    stamps = strsplit(name, '__');
    names = cell(size(stamps));
    for i = 1:numel(stamps)
        names{i} = stampedName(stamps{i});
        if isempty(names{i})
            return;
        end
    end
    name = strjoin(names, '+');
    vouched = {name};
end

function name = stampedName(stamp)
    % The stage kind that `stamp` is a generated stamp of, or ''. A stamp is the name of its kind;
    % kinds with stamped options append an underscore and the option values.
    kinds = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
    bare = {FeatureName.PLAIN.value, FeatureName.PERMUTED.value, FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value, FeatureName.CUSTOM.value};
    name = '';
    for i = 1:numel(kinds)
        if strcmp(stamp, kinds{i}) || (~ismember(kinds{i}, bare) && startsWith(stamp, [kinds{i}, '_']))
            name = kinds{i};
            return;
        end
    end
end
