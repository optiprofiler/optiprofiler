classdef FeatureName
%FEATURENAME enumerates all the possible feature names.
    
    enumeration
        PLAIN ('plain')
        PERTURBED_X0 ('perturbed_x0')
        NOISY ('noisy')
        TRUNCATED ('truncated')
        PERMUTED ('permuted')
        ROTATED ('rotated')
        BADLY_SCALED ('badly_scaled')
        RANDOM_NAN ('random_nan')
        UNRELAXABLE_CONSTRAINTS ('unrelaxable_constraints')
        CUSTOM ('custom')
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