classdef FeatureOptionKey
%FEATUREOPTIONKEY enumerates options for defining features.
    
    enumeration
        N_RUNS ('n_runs')
        DISTRIBUTION ('distribution')
        NOISE_LEVEL ('noise_level')
        NOISE_TYPE ('noise_type')
        SIGNIFICANT_DIGITS ('significant_digits')
        PERTURBED_TRAILING_ZEROS ('perturbed_trailing_zeros')
        ROTATED ('rotated')
        CONDITION_FACTOR ('condition_factor')
        RATE_NAN ('rate_nan')
        UNRELAXABLE_BOUNDS ('unrelaxable_bounds')
        UNRELAXABLE_LINEAR_CONSTRAINTS ('unrelaxable_linear_constraints')
        UNRELAXABLE_NONLINEAR_CONSTRAINTS ('unrelaxable_nonlinear_constraints')
        MESH_SIZE ('mesh_size')
        IS_TRUTH ('is_truth')
        MOD_X0 ('mod_x0')
        MOD_AFFINE ('mod_affine')
        MOD_BOUNDS ('mod_bounds')
        MOD_LINEAR_UB ('mod_linear_ub')
        MOD_LINEAR_EQ ('mod_linear_eq')
        MOD_FUN ('mod_fun')
        MOD_CUB ('mod_cub')
        MOD_CEQ ('mod_ceq')
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