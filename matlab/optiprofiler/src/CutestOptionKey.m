classdef CutestOptionKey
%CUTESTOPTIONKEY enumerates options for selecting problems and constructing the problem set.
    
    enumeration
        P_TYPE ('p_type')
        MINDIM ('mindim')
        MAXDIM ('maxdim')
        MINCON ('mincon')
        MAXCON ('maxcon')
        EXCLUDELIST ('excludelist')
    end
    properties
        value
    end
    methods
        function obj = CutestOptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end