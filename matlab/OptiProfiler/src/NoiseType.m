classdef NoiseType
%NOISETYPE emumerates the types of noise that can be added to the objective function.
    
    enumeration
        ABSOLUTE ('absolute')
        RELATIVE ('relative')
    end
    properties
        value
    end
    methods
        function obj = NoiseType(inputValue)
            obj.value = inputValue;
        end
    end
end