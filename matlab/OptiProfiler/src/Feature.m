classdef Feature < handle
    %FEATURE used to modify the objective function.

    properties (GetAccess = public, SetAccess = private)

        name
        options
        
    end

    methods

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
        function obj = Feature(varargin)
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
                if ~ismember(optionKeys{i}, validOptionKeys)
                    error("Unknown option: " + optionKeys{i} + ".")
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
                for i = 1:numel(optionKeys)
                    if ~ismember(optionKeys{i}, known_options)
                        error("Option " + optionKeys{i} + " is not valid for feature " + obj.name + ".")
                    end
                end





            end
        end

    end
end