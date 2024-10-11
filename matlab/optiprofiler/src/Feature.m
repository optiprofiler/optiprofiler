classdef Feature < handle
%FEATURE is a class that represents different kinds of modifications to the optimization problem.

    properties (GetAccess = public, SetAccess = private)

        name
        options
        
    end

    methods

        function obj = Feature(feature_name, varargin)
            %{
            Initialize a feature.

            Parameters
            ----------
            feature_name : str
                Name of the feature.

            Other Parameters
            ----------------
            distribution : callable, optional
                Distribution used in the noisy and perturbed_x0 features.
            modifier : callable, optional
                Modifier used in the custom feature.
            n_runs : int, optional
                Number of runs of the feature.
            rate_nan : int or float, optional
                Rate of NaNs of the random_nan feature.
            significant_digits : int, optional
                Number of significant digits of the truncated feature.
            noise_type : str, optional
                Type of the noisy and perturbed_x0 features.
            %}

            if ~ischarstr(feature_name)
                error("MATLAB:Feature:FeaturenameNotString", "The feature name must be a char or a string.")
            end
            % Convert the feature name to lowercase characters.
            obj.name = char(lower(feature_name));

            obj.options = struct();
            if nargin > 1
                if isstruct(varargin{1})
                    obj.options = varargin{1};
                elseif mod(length(varargin), 2) ~= 0
                    error("MATLAB:Feature:InvalidNumberOfArguments", "The input of the feature options must be in pairs.")
                else
                    for i = 1:2:length(varargin)
                        obj.options.(lower(varargin{i})) = varargin{i+1};
                    end
                end
            end
    
            % Check whether the feature is valid.
            % validFeatureNames = {enumeration('FeatureName').value};  % Only for MATLAB R2021b and later
            validFeatureNames = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
            if ~ismember(obj.name, validFeatureNames)
                error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end
            optionKeys = fieldnames(obj.options);
            validOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('FeatureOptionKey')), 'UniformOutput', false);
            for i = 1:numel(optionKeys)
                key = optionKeys{i};
                if ~ismember(key, validOptionKeys)
                    error("MATLAB:Feature:UnknownOption", "Unknown option: " + key + ".")
                end

                % Check whether the option is valid for the feature.
                known_options = {FeatureOptionKey.N_RUNS.value};
                switch obj.name
                    case FeatureName.CUSTOM.value
                        known_options = [known_options, {FeatureOptionKey.MODIFIER.value}];
                    case FeatureName.NOISY.value
                        known_options = [known_options, {FeatureOptionKey.DISTRIBUTION.value, FeatureOptionKey.NOISE_LEVEL.value, FeatureOptionKey.NOISE_TYPE.value}];
                    case FeatureName.PERTURBED_X0.value
                        known_options = [known_options, {FeatureOptionKey.DISTRIBUTION.value, FeatureOptionKey.NOISE_LEVEL.value}];
                    case FeatureName.RANDOM_NAN.value
                        known_options = [known_options, {FeatureOptionKey.RATE_NAN.value}];
                    case FeatureName.TRUNCATED.value
                        known_options = [known_options, {FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value, FeatureOptionKey.SIGNIFICANT_DIGITS.value}];
                    case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                        known_options = [known_options, {FeatureOptionKey.UNRELAXABLE_BOUNDS.value, FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value, FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value}];
                    case FeatureName.LINEARLY_TRANSFORMED.value
                        known_options = [known_options, {FeatureOptionKey.ROTATED.value, FeatureOptionKey.INVERTIBLE_TRANSFORMATION.value, FeatureOptionKey.CONDITION_NUMBER.value}];
                    case FeatureName.PERMUTED.value
                        % Do nothing
                    case FeatureName.PLAIN.value
                        % Do nothing
                    otherwise
                        error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
                end
                if ~ismember(key, known_options)
                    error("MATLAB:Feature:InvalidOptionForFeature", "Option " + key + " is not valid for feature " + obj.name + ".")
                end

                % Check whether the option type is valid.
                switch key
                    case FeatureOptionKey.N_RUNS.value
                        if ~isintegerscalar(obj.options.(key)) || obj.options.(key) <= 0
                            error("MATLAB:Feature:n_runs_NotPositiveInteger", "Option " + key + " must be a positive integer.")
                        end
                    case FeatureOptionKey.MODIFIER.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:modifier_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.DISTRIBUTION.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:distribution_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.RATE_NAN.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) < 0.0 || obj.options.(key) > 1.0
                            error("MATLAB:Feature:rate_nan_NotBetween_0_1", "Option " + key + " must be a real number between 0 and 1.")
                        end
                    case FeatureOptionKey.SIGNIFICANT_DIGITS.value
                        if ~isintegerscalar(obj.options.(key)) || obj.options.(key) <= 0
                            error("MATLAB:Feature:significant_digits_NotPositiveInteger", "Option " + key + " must be a positive integer.")
                        end
                    case FeatureOptionKey.NOISE_LEVEL.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) <= 0.0
                            error("MATLAB:Feature:noise_level_NotPositive", "Option " + key + " must be a positive real number.")
                        end
                    case FeatureOptionKey.NOISE_TYPE.value
                        validNoiseTypes = cellfun(@(x) x.value, num2cell(enumeration('NoiseType')), 'UniformOutput', false);
                        if ~ischarstr(obj.options.(key)) || ~ismember(obj.options.(key), validNoiseTypes)
                            error("MATLAB:Feature:noise_type_InvalidInput", "Option " + key + " must be '" + NoiseType.ABSOLUTE.value + "' or '" + NoiseType.RELATIVE.value + "' or '" + NoiseType.MIXED.value + "'.")
                        end
                    case FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:perturbed_trailing_zeros_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.ROTATED.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:rotated_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.INVERTIBLE_TRANSFORMATION.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:invertible_transformation_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.CONDITION_NUMBER.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:condition_number_InvalidInput", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_BOUNDS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_bounds_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_linear_constraints_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_nonlinear_constraints_NotLogical", "Option " + key + " must be a logical.")
                        end
                end
            end

            % Set default options.
            obj = set_default_options(obj);
        end

        function is_stochastic = isStochastic(obj)
            stochasticFeatures = {FeatureName.CUSTOM.value, FeatureName.NOISY.value, FeatureName.PERMUTED.value, FeatureName.LINEARLY_TRANSFORMED.value, FeatureName.PERTURBED_X0.value, FeatureName.RANDOM_NAN.value, FeatureName.TRUNCATED.value};
            is_stochastic = ismember(obj.name, stochasticFeatures);
        end

        function f = modifier(obj, x, f, seed, maxcv_bounds, maxcv_linear, maxcv_nonlinear)
            %{
            Modify the objective function value.

            Parameters
            ----------
            x : `numpy.ndarray`, shape (n,)
                Point at which the objective function is evaluated.
            f : float
                Objective function value at `x`.
            seed : int, optional
                Seed used to generate random numbers.

            Returns
            -------
            float
                Modified objective function value.
            %}

            if nargin < 4
                seed = 'shuffle';
            end
            if nargin < 5
                maxcv_bounds = 0;
            end
            if nargin < 6
                maxcv_linear = 0;
            end
            if nargin < 7
                maxcv_nonlinear = 0;
            end
            

            % Convert x into a cell array
            xCell = num2cell(x);

            switch obj.name
                case FeatureName.CUSTOM.value
                    f = obj.options.(FeatureOptionKey.MODIFIER.value)(x, f, seed);
                case FeatureName.NOISY.value
                    rand_stream = obj.default_rng(seed, f, obj.options.(FeatureOptionKey.NOISE_LEVEL.value), sum(double(obj.options.(FeatureOptionKey.NOISE_TYPE.value))), xCell{:});
                    if obj.options.(FeatureOptionKey.NOISE_TYPE.value) == NoiseType.ABSOLUTE.value
                        f = f + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream, 1);
                    elseif obj.options.(FeatureOptionKey.NOISE_TYPE.value) == NoiseType.RELATIVE.value
                        f = f * (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream, 1));
                    else
                        % We need the distribution of the noise to be symmetric with respect to 0.
                        f = f + max(1, abs(f)) * obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream, 1);
                    end
                case FeatureName.RANDOM_NAN.value
                    rand_stream = obj.default_rng(seed, f, obj.options.(FeatureOptionKey.RATE_NAN.value), xCell{:});
                    if rand_stream.rand() < obj.options.(FeatureOptionKey.RATE_NAN.value)
                        f = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    rand_stream = obj.default_rng(seed, f, obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value), xCell{:});
                    if f == 0
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(f))) - 1;
                    end
                    digits = double(digits);
                    % Round f to the desired number of significant digits. (We can also use 'fix(x)' to round towards zero.)
                    f = round(f, digits);
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                        if f >= 0
                            f = f + rand_stream.rand() * 10 ^ (-digits);
                        else
                            f = f - rand_stream.rand() * 10 ^ (-digits);
                        end
                    end
                case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                    if obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) && maxcv_bounds > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value) && maxcv_linear > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value) && maxcv_nonlinear > 0
                        f = Inf;
                    end
                case FeatureName.PERMUTED.value
                    % Do nothing
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Do nothing
                case FeatureName.PLAIN.value
                    % Do nothing
                case FeatureName.PERTURBED_X0.value
                    % Do nothing
                otherwise
                    error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end
        end

        function obj = set_default_options(obj)
            % Set default options.
            switch obj.name
                case FeatureName.PLAIN.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.CUSTOM.value
                    if ~isfield(obj.options, FeatureOptionKey.MODIFIER.value)
                        error("MATLAB:Feature:MissingModifier", "When using a custom feature, you must specify the " + FeatureOptionKey.MODIFIER.value + " option.");
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.NOISY.value
                    if ~isfield(obj.options, FeatureOptionKey.DISTRIBUTION.value)
                        obj.options.(FeatureOptionKey.DISTRIBUTION.value) = @obj.default_distribution;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_LEVEL.value)
                        obj.options.(FeatureOptionKey.NOISE_LEVEL.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_TYPE.value)
                        obj.options.(FeatureOptionKey.NOISE_TYPE.value) = NoiseType.MIXED.value;
                    end
                case FeatureName.PERMUTED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                case FeatureName.LINEARLY_TRANSFORMED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.ROTATED.value)
                        obj.options.(FeatureOptionKey.ROTATED.value) = true;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.INVERTIBLE_TRANSFORMATION.value)
                        obj.options.(FeatureOptionKey.INVERTIBLE_TRANSFORMATION.value) = @obj.default_invertible_transformation;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.CONDITION_NUMBER.value)
                        obj.options.(FeatureOptionKey.CONDITION_NUMBER.value) = @(n) 1;
                    end
                case FeatureName.PERTURBED_X0.value
                    if ~isfield(obj.options, FeatureOptionKey.DISTRIBUTION.value)
                        obj.options.(FeatureOptionKey.DISTRIBUTION.value) = @obj.default_distribution;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_LEVEL.value)
                        obj.options.(FeatureOptionKey.NOISE_LEVEL.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                case FeatureName.RANDOM_NAN.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.RATE_NAN.value)
                        obj.options.(FeatureOptionKey.RATE_NAN.value) = 0.05;
                    end
                case FeatureName.TRUNCATED.value
                    if ~isfield(obj.options, FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                        obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value) = true;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                            obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                        else
                            obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                        end
                    end
                    if ~isfield(obj.options, FeatureOptionKey.SIGNIFICANT_DIGITS.value)
                        obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) = 6;
                    end
                case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.UNRELAXABLE_BOUNDS.value)
                        obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) = false;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value)
                        obj.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value) = false;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value)
                        obj.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value) = false;
                    end
                otherwise
                    error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end
        end

    end

    methods (Static)
        function value = default_distribution(rand_stream, n)
            if nargin < 2
                n = 1;
            end
            value = rand_stream.randn(n, 1);
        end

        function [transform, inverse] = default_invertible_transformation(rand_stream, n)
            % Define the default invertible transformation for the 'linearly_transformed' feature.

            % We generate a random orthogonal matrix Q following the uniform distribution on SO(n).
            % The method refers to a note written by the late Professor Nicholas Higham, see:
            % https://nhigham.com/2020/04/22/what-is-a-random-orthogonal-matrix/
            % and an answer from MathStackExchange, see:
            % https://math.stackexchange.com/a/4891933/1088047

            [Q, R] = qr(rand_stream.randn(n));
            Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
            transform = Q;
            inverse = Q';
        end

        function rand_stream = default_rng(seed, varargin)
            % Generate a random number generator.
            %
            % Parameters
            % ----------
            % seed : double, but default one is 'shuffle'
            %     Seed used to generate an initial random number generator.
            % varargin : array of double
            %     Arguments used to generate the returned random number generator.
            %
            % Returns
            % -------
            % rand_stream : RandStream
            %     Random number generator.

            % Create an initial rand_stream with the given seed
            if ~strcmp(seed, 'shuffle')
                if isnan(seed)
                    seed = 0;
                end
                if ~isrealscalar(seed)
                    error("MATLAB:Feature:SeedNotEvenReal", "The input for seed should at least be a real number.")
                end
                seed = mod(floor(seed), 2^32);
            end
            rand_stream = RandStream('mt19937ar', 'Seed', seed);

            % Convert all elements in varargin to double
            varargin = cellfun(@double, varargin, 'UniformOutput', false);

            % Generate a new seed based on the initial rand_stream and the additional arguments
            newSeed = abs(sin(1e5 * randn(rand_stream, 1)) + sum(sin(1e5 * prod(cellfun(@(x) x, varargin))))) * 1e9;
            newSeed = mod(floor(newSeed), 2^32);

            % Create a new rand_stream with the new seed
            if isnan(newSeed)
                newSeed = 0;
            end
            rand_stream = RandStream('mt19937ar', 'Seed', floor(newSeed));
        end
    end
end