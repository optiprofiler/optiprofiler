classdef OtherOptionKey
%OTHEROPTIONKEY enumerates all other options except for Cutest, features, or profiles.
    
    enumeration
        PROBLEM ('problem')
        CUTEST_PROBLEM_NAMES ('cutest_problem_names')
        CUSTOM_PROBLEM_LOADER ('custom_problem_loader')
        CUSTOM_PROBLEM_NAMES ('custom_problem_names')
    end
    properties
        value
    end
    methods
        function obj = OtherOptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end