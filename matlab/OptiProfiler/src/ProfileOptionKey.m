classdef ProfileOptionKey
    %PROFILEOPTIONKEY of the available options.
    %
    
    enumeration
        N_JOBS ('n_jobs')
    end
    properties
        value
    end
    methods
        function obj = ProfileOptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end