function name = archivedFeatureName(result)
%ARCHIVEDFEATURENAME returns the name of the feature of an archived library record, for display.
%   A saved experiment has the `feature_stamp` of its feature. A generated stamp is the feature name,
%   followed by an underscore and option values for the features that have options in their stamps
%   ('linearly_transformed_rotated'), so the name is read from it. OptiProfiler 2 joins the stamps of
%   composed features with a double underscore; their names are joined with '+', as there. A stamp of
%   another form (it can be set by the user) is returned as it is, since it is the only name that a
%   saved experiment has. The name is empty when the record has no stamp.

    name = '';
    if ~isfield(result, 'feature_stamp') || ~ischarstr(result.feature_stamp) || ~isscalar(string(result.feature_stamp)) ...
            || strlength(result.feature_stamp) == 0
        return;
    end
    name = char(result.feature_stamp);
    kinds = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
    bare = {FeatureName.PLAIN.value, FeatureName.PERMUTED.value, FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value, FeatureName.CUSTOM.value};
    stamps = strsplit(name, '__');
    names = cell(size(stamps));
    for i = 1:numel(stamps)
        for j = 1:numel(kinds)
            if strcmp(stamps{i}, kinds{j}) || (~ismember(kinds{j}, bare) && startsWith(stamps{i}, [kinds{j}, '_']))
                names{i} = kinds{j};
                break;
            end
        end
        if isempty(names{i})
            return;
        end
    end
    name = strjoin(names, '+');
end
