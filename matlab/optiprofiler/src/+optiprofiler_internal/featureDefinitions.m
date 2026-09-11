function definition = featureDefinitions(name, options)
%FEATUREDEFINITIONS Pure local schemas and literal experiment replicate hints.
% The optional OPTIONS input is already normalized; it only resolves the hint
% and stochastic predicate. This function never constructs runtime objects.
    names = {'plain','perturbed_x0','noisy','truncated','permuted','linearly_transformed', ...
        'random_nan','unrelaxable_constraints','nonquantifiable_constraints','quantized','custom'};
    if nargin == 0
        definition = repmat(optiprofiler_internal.featureDefinitions('plain'),1,numel(names));
        for k = 1:numel(names), definition(k) = optiprofiler_internal.featureDefinitions(names{k}); end
        return;
    end
    definition = struct('name', name, 'code', [], 'local_keys', {{}}, ...
        'defaults', struct(), 'replicate_hint', 1, 'is_stochastic', false);
    switch name
        case 'plain'
        case 'perturbed_x0'
            definition.code = 1;
            definition.local_keys = {'distribution','perturbation_level'};
            definition.defaults = struct('distribution','spherical','perturbation_level',1e-3);
            definition.replicate_hint = 5; definition.is_stochastic = true;
        case 'noisy'
            definition.code = 2;
            definition.local_keys = {'distribution','noise_level','noise_type','noise_mode','noise_map'};
            definition.defaults = struct('noise_mode','random','distribution','gaussian', ...
                'noise_map','chebyshev','noise_level',1e-3,'noise_type','mixed');
        case 'truncated'
            definition.code = 3;
            definition.local_keys = {'perturbed_trailing_digits','significant_digits'};
            definition.defaults = struct('perturbed_trailing_digits',false,'significant_digits',6);
        case 'permuted'
            definition.code = 4;
            definition.replicate_hint = 5; definition.is_stochastic = true;
        case 'linearly_transformed'
            definition.code = 5;
            definition.local_keys = {'rotated','condition_factor'};
            definition.defaults = struct('rotated',true,'condition_factor',0);
            definition.replicate_hint = 5;
        case 'random_nan'
            definition.code = 6;
            definition.local_keys = {'nan_rate'};
            definition.defaults = struct('nan_rate',.05);
            definition.replicate_hint = 5; definition.is_stochastic = true;
        case 'unrelaxable_constraints'
            definition.code = 7;
            definition.local_keys = {'unrelaxable_bounds','unrelaxable_linear_constraints','unrelaxable_nonlinear_constraints'};
            definition.defaults = struct('unrelaxable_bounds',true, ...
                'unrelaxable_linear_constraints',false,'unrelaxable_nonlinear_constraints',false);
        case 'nonquantifiable_constraints'
            definition.code = 8;
        case 'quantized'
            definition.code = 9;
            definition.local_keys = {'mesh_size','mesh_type','ground_truth'};
            definition.defaults = struct('mesh_size',1e-3,'mesh_type','absolute','ground_truth',true);
        case 'custom'
            definition.code = 10;
            definition.local_keys = {'mod_x0','mod_bounds','mod_linear_ub','mod_linear_eq', ...
                'mod_affine','mod_fun','mod_cub','mod_ceq'};
            definition.is_stochastic = true;
        otherwise
            error('MATLAB:Feature:UnknownFeature', 'Unknown feature: %s.', name);
    end
    if nargin < 2, options = definition.defaults; end
    if strcmp(name,'noisy')
        definition.is_stochastic = strcmp(options.noise_mode,'random');
        if definition.is_stochastic, definition.replicate_hint = 5; end
    elseif strcmp(name,'truncated')
        definition.is_stochastic = options.perturbed_trailing_digits;
        if definition.is_stochastic, definition.replicate_hint = 5; end
    elseif strcmp(name,'linearly_transformed')
        definition.is_stochastic = options.rotated;
    end
end
