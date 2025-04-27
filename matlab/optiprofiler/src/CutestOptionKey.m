classdef CutestOptionKey
%CUTESTOPTIONKEY enumerates options for selecting problems and constructing the problem set.
    
    enumeration
        PTYPE ('ptype')
        MINDIM ('mindim')
        MAXDIM ('maxdim')
        MINB ('minb')
        MAXB ('maxb')
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