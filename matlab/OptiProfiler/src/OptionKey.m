classdef OptionKey
    %OPTIONKEY of the available options.
    %
    
    enumeration
        DISTRIBUTION ('distribution')
        MODIFIER ('modifier')
        N_RUNS ('n_runs')
        ORDER ('order')
        PARAMETER ('parameter')
        RATE_ERROR ('rate_error')
        RATE_NAN ('rate_nan')
        SIGNIFICANT_DIGITS ('significant_digits')
        TYPE ('type')
    end
    properties
        value
    end
    methods
        function obj = OptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end