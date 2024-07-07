classdef FeatureOptionKey
%FEATUREOPTIONKEY enumerates options for defining features.
    
    enumeration
        N_RUNS ('n_runs')
        DISTRIBUTION ('distribution')
        NOISE_LEVEL ('noise_level')
        NOISE_TYPE ('noise_type')
        SIGNIFICANT_DIGITS ('significant_digits')
        UNRELAXABLE_BOUNDS ('unrelaxable_bounds')
        UNRELAXABLE_LINEAR_CONSTRAINTS ('unrelaxable_linear_constraints')
        UNRELAXABLE_NONLINEAR_CONSTRAINTS ('unrelaxable_nonlinear_constraints')
        PERTURBED_TRAILING_ZEROS ('perturbed_trailing_zeros')
        RATE_NAN ('rate_nan')
        MODIFIER ('modifier')
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