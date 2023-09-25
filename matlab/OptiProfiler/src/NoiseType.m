classdef NoiseType
    %NOISETYPE of the available noise types.

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