classdef FeatureName
%FEATURENAME enumerates all the possible feature names.
    
    enumeration
        CUSTOM ('custom')
        NOISY ('noisy')
        PERMUTED ('permuted')
        PERTURBED_X0 ('perturbed_x0')
        PLAIN ('plain')
        RANDOM_NAN ('random_nan')
        TRUNCATED ('truncated')
        UNRELAXABLE_CONSTRAINTS ('unrelaxable_constraints')
    end
    properties
        value
    end
    methods
        function obj = FeatureName(inputValue)
            obj.value = inputValue;
        end
    end
end