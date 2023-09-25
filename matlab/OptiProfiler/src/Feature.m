classdef Feature < handle
    %FEATURE used to modify the objective function.
    %
    properties (GetAccess = public, SetAccess = private)

        name
        options
        
    end

    methods

        function obj = Feature(varargin)
            %{
            Initialize a feature.

            Parameters
            ----------
            feature_name : str
                Name of the feature.

            Other Parameters
            ----------------
            distribution : callable, optional
                Distribution used in the noisy feature.
            modifier : callable, optional
                Modifier used in the custom feature.
            order : int or float, optional
                Order of the regularized feature.
            parameter : int or float, optional
                Regularization parameter of the regularized feature.
            rate_error : int or float, optional
                Rate of errors of the tough feature.
            rate_nan : int or float, optional
                Rate of NaNs of the tough feature.
            significant_digits : int, optional
                Number of significant digits of the truncated feature.
            type : str, optional
                Type of the noisy feature.
            %}
            if isempty(varargin)
                error("Feature name must be provided.")
            end
            obj.name = lower(varargin{1});
            obj.options = struct();
            for i = 2:2:length(varargin)
                obj.options.(lower(varargin{i})) = varargin{i+1};
            end
    
            % Check whether the feature is valid.
            validFeatureNames = {enumeration('FeatureName').value};
            if ~ismember(obj.name, validFeatureNames)
                error("Unknown feature: " + obj.name + ".")
            end
            optionKeys = fieldnames(obj.options);
            validOptionKeys = {enumeration('OptionKey').value};
            for i = 1:numel(optionKeys)
                key = optionKeys{i};
                if ~ismember(key, validOptionKeys)
                    error("Unknown option: " + key + ".")
                end

                % Check whether the option is valid for the feature.
                known_options = {OptionKey.N_RUNS.value};
                switch obj.name
                    case FeatureName.CUSTOM.value
                        known_options = [known_options, {OptionKey.MODIFIER.value}];
                    case FeatureName.NOISY.value
                        known_options = [known_options, {OptionKey.DISTRIBUTION.value, OptionKey.TYPE.value}];
                    case FeatureName.REGULARIZED.value
                        known_options = [known_options, {OptionKey.ORDER.value, OptionKey.PARAMETER.value}];
                    case FeatureName.TOUGH.value
                        known_options = [known_options, {OptionKey.RATE_ERROR.value, OptionKey.RATE_NAN.value}];
                    case FeatureName.TRUNCATED.value
                        known_options = [known_options, {OptionKey.SIGNIFICANT_DIGITS.value}];
                    case FeatureName.PLAIN.value
                        % Do nothing
                    otherwise
                        error("Unknown feature: " + obj.name + ".")
                end
                if ~ismember(key, known_options)
                    error("Option " + key + " is not valid for feature " + obj.name + ".")
                end

                % Check whether the option type is valid.
                switch key
                    case OptionKey.N_RUNS.value
                        if isfloat(obj.options.(key)) && (obj.options.(key) == round(obj.options.(key)))
                            obj.options.(key) = int32(obj.options.(key));
                        end
                        if ~isinteger(obj.options.(key)) || obj.options.(key) <= 0
                            error("Option " + key + " must be a positive integer.")
                        end
                    case OptionKey.MODIFIER.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("Option " + key + " must be a function handle.")
                        end
                    case OptionKey.DISTRIBUTION.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("Option " + key + " must be a function handle.")
                        end
                    case OptionKey.ORDER.value
                        if ~isnumic(obj.options.(key))
                            error("Option " + key + " must be a number.")
                        end
                    case OptionKey.PARAMETER.value
                        if ~isnumic(obj.options.(key)) || obj.options.(key) < 0.0
                            error("Option " + key + " must be a nonnegative number.")
                        end
                    case {OptionKey.RATE_ERROR.value, OptionKey.RATE_NAN.value}
                        if ~isnumic(obj.options.(key)) || obj.options.(key) < 0.0 || obj.options.(key) > 1.0
                            error("Option " + key + " must be a number between 0 and 1.")
                        end
                    case OptionKey.SIGNIFICANT_DIGITS.value
                        if isfloat(obj.options.(key)) && (obj.options.(key) == round(obj.options.(key)))
                            obj.options.(key) = int32(obj.options.(key));
                        end
                        if ~isinteger(obj.options.(key)) || obj.options.(key) <= 0
                            error("Option " + key + " must be a positive integer.")
                        end
                    case OptionKey.TYPE.value
                        validNoiseTypes = {enumeration('NoiseType').value};
                        if ~(ischar(obj.options.(key)) || isstring(obj.options.(key))) || ~ismember(obj.options.(key), validNoiseTypes)
                            error("Option " + key + " must be either '" + NoiseType.ABSOLUTE.value + "' or '" + NoiseType.RELATIVE.value + "'.")
                        end
                end
            end

            % Set default options.
            obj = set_default_options(obj);
        end

        function f = modifier(obj, x, f, seed)
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
                seed = [];
            end

            % Convert x into a cell array
            xCell = num2cell(x);

            switch obj.name
                case FeatureName.PLAIN.value
                    % Do nothing
                case FeatureName.CUSTOM.value
                    f = obj.options.(OptionKey.MODIFIER.value)(x, f, seed);
                case FeatureName.NOISY.value
                    rng = obj.default_rng(seed, f, sum(double(obj.options.(OptionKey.TYPE.value))), xCell{:});
                    if obj.options.(OptionKey.TYPE.value) == NoiseType.ABSOLUTE.value
                        f = f + obj.options.(OptionKey.DISTRIBUTION.value)(rng);
                    else
                        f = f * (1.0 + obj.options.(OptionKey.DISTRIBUTION.value)(rng));
                    end
                case FeatureName.REGULARIZED.value
                    f = f + obj.options.(OptionKey.PARAMETER.value) * norm(x, obj.options.(OptionKey.ORDER.value));
                case FeatureName.TOUGH.value
                    rng = obj.default_rng(seed, f, obj.options.(OptionKey.RATE_ERROR.value), obj.options.(OptionKey.RATE_NAN.value), xCell{:});
                    if rng.rand() < obj.options.(OptionKey.RATE_ERROR.value)
                        error("Runtime error.")
                    elseif rng.rand() < obj.options.(OptionKey.RATE_NAN.value)
                        f = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    rng = obj.default_rng(seed, f, obj.options.(OptionKey.SIGNIFICANT_DIGITS.value), xCell{:});
                    if f == 0.0
                        digits = obj.options.(OptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(OptionKey.SIGNIFICANT_DIGITS.value) - int32(floor(log10(abs(f)))) - 1;
                    end
                    if f >= 0.0
                        f = round(f * 10^digits) / 10^digits + rng.rand() * 10^(-digits);
                    else
                        f = round(f * 10^digits) / 10^digits - rng.rand() * 10^(-digits);
                    end
                otherwise
                    error("Unknown feature: " + obj.name + ".")
            end
        end

        function obj = set_default_options(obj)
            % Set default options.
            switch obj.name
                case FeatureName.PLAIN.value
                    if ~isfield(obj.options, OptionKey.N_RUNS.value)
                        obj.options.(OptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.CUSTOM.value
                    if ~isfield(obj.options, OptionKey.MODIFIER.value)
                        error("When using a custom feature, you must specify the " + OptionKey.MODIFIER.value + " option.");
                    end
                    if ~isfield(obj.options, 'n_runs')
                        obj.options.(OptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.NOISY.value
                    if ~isfield(obj.options, OptionKey.DISTRIBUTION.value)
                        obj.options.(OptionKey.DISTRIBUTION.value) = @(rng) 1e-3 * rng.randn();
                    end
                    if ~isfield(obj.options, OptionKey.N_RUNS.value)
                        obj.options.(OptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, OptionKey.TYPE.value)
                        obj.options.(OptionKey.TYPE.value) = NoiseType.RELATIVE.value;
                    end
                case FeatureName.REGULARIZED.value
                    if ~isfield(obj.options, OptionKey.N_RUNS.value)
                        obj.options.(OptionKey.N_RUNS.value) = 1;
                    end
                    if ~isfield(obj.options, OptionKey.ORDER.value)
                        obj.options.(OptionKey.ORDER.value) = 2;
                    end
                    if ~isfield(obj.options, OptionKey.PARAMETER.value)
                        obj.options.(OptionKey.PARAMETER.value) = 1.0;
                    end
                case FeatureName.TOUGH.value
                    if ~isfield(obj.options, OptionKey.N_RUNS.value)
                        obj.options.(OptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, OptionKey.RATE_ERROR.value)
                        obj.options.(OptionKey.RATE_ERROR.value) = 0.0;
                    end
                    if ~isfield(obj.options, OptionKey.RATE_NAN.value)
                        obj.options.(OptionKey.RATE_NAN.value) = 0.05;
                    end
                case FeatureName.TRUNCATED.value
                    if ~isfield(obj.options, OptionKey.N_RUNS.value)
                        obj.options.(OptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, OptionKey.SIGNIFICANT_DIGITS.value)
                        obj.options.(OptionKey.SIGNIFICANT_DIGITS.value) = 6;
                    end
                otherwise
                    error("Unknown feature: " + obj.name + ".")
            end
        end

    end

    methods (Static)
        function rng = default_rng(seed, varargin)
            % Generate a random number generator.
            %
            % Parameters
            % ----------
            % seed : double
            %     Seed used to generate an initial random number generator.
            % varargin : array of double
            %     Arguments used to generate the returned random number generator.
            %
            % Returns
            % -------
            % rng : RandStream
            %     Random number generator.

            % Create an initial RNG with the given seed
            rng = RandStream('mt19937ar', 'Seed', seed);

            % Generate a new seed based on the initial RNG and the additional arguments
            newSeed = abs(sin(1e5 * randn(rng, 1)) + sum(sin(1e5 * varargin{:}))) * 1e9;

            % Create a new RNG with the new seed
            rng = RandStream('mt19937ar', 'Seed', floor(newSeed));
        end
    end
end