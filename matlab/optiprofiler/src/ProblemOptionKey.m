classdef ProblemOptionKey
%CUTESTOPTIONKEY enumerates options for selecting problems and constructing the problem set.

    enumeration
        PLIBS ('plibs')
        PTYPE ('ptype')
        MINDIM ('mindim')
        MAXDIM ('maxdim')
        MINB ('minb')
        MAXB ('maxb')
        MINCON ('mincon')
        MAXCON ('maxcon')
        EXCLUDELIST ('excludelist')
        PROBLEM_NAMES ('problem_names')
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