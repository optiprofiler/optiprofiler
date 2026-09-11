classdef LegacyFeatureEnvelope
%LEGACYFEATUREENVELOPE Retained old native data, deliberately not a Feature.
% Native load is trusted input. Convert through importLegacyFeature before
% execution; the old resolved repetition count does not belong to new Feature.
    properties (SetAccess = private)
        saved_state
    end
    methods
        function obj = LegacyFeatureEnvelope(saved)
            if ~isa(saved,'struct') || ~isscalar(saved) || ...
                    ~isequal(sort(fieldnames(saved)),sort({'name';'options'}))
                error('MATLAB:Feature:InvalidNativeState','Unsupported historical Feature layout.');
            end
            obj.saved_state = saved;
        end
    end
end
