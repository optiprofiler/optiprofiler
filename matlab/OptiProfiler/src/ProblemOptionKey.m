classdef ProblemOptionKey
    %PROBLEMOPTIONKEY Problem options for CUTEst problems
    %

    enumeration
        N_MIN ('mindim')
        N_MAX ('maxdim')
        M_MIN ('mincon')
        M_MAX ('maxcon')
    end
    properties
        value
    end
    methods
        function obj = ProblemOptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end