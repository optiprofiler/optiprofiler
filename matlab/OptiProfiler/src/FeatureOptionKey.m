classdef FeatureOptionKey
%FEATUREOPTIONKEY enumerates options for defining features.
    
    enumeration
        DISTRIBUTION ('distribution')
        MODIFIER ('modifier')
        N_RUNS ('n_runs')
        RATE_NAN ('rate_nan')
        PERTURBED_TRAILING_ZEROS ('perturbed_trailing_zeros')
        SIGNIFICANT_DIGITS ('significant_digits')
        TYPE ('type')
        UNRELAXABLE_BOUNDS ('unrelaxable_bounds')
        UNRELAXABLE_LINEAR_CONSTRAINTS ('unrelaxable_linear_constraints')
        UNRELAXABLE_NONLINEAR_CONSTRAINTS ('unrelaxable_nonlinear_constraints')
    end
    properties
        value
    end
    methods
        function obj = FeatureOptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end