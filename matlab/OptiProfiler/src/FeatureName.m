classdef FeatureName
%FEATURENAME enumerates all the possible feature names.
    
    enumeration
        CUSTOM ('custom')
        NOISY ('noisy')
        PLAIN ('plain')
        RANDOMIZE_X0 ('randomize_x0')
        REGULARIZED ('regularized')
        TOUGH ('tough')
        TRUNCATED ('truncated')
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