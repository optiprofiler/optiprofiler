classdef Feature < handle
%FEATURE Reusable, locally configured ordered problem transformations.
% STAGES and DECLARED return value copies; user function handles remain native
% references and are not invoked during normalization. Experiment repetitions
% belong to benchmark options, never to this specification.
    properties (Access = private)
        specification_
    end
    properties (Dependent, SetAccess = private)
        stages
        declared
        name
        declared_name
        is_stochastic
        is_identity
        specification_version
        options
    end
    methods
        function obj = Feature(input, varargin)
            if nargin == 0
                error('MATLAB:Feature:MissingSpecification','Feature requires an explicit specification, such as Feature(''plain'').');
            end
            if strcmp(class(input), 'Feature') && isscalar(input)
                if ~isempty(varargin)
                    error('MATLAB:Feature:StructuredOverrides', ...
                        'Feature-object input accepts no extra options; n_runs belongs in benchmark options.');
                end
                obj = input;
                return;
            end
            obj.specification_ = optiprofiler_internal.normalizeFeatureSpecification(input, varargin{:});
        end
        function value = get.stages(obj), value = obj.specification_.stages; end
        function value = get.declared(obj), value = obj.specification_.declared; end
        function value = get.name(obj), value = obj.specification_.name; end
        function value = get.declared_name(obj), value = obj.specification_.declared_name; end
        function value = get.is_stochastic(obj), value = obj.specification_.is_stochastic; end
        function value = get.is_identity(obj), value = obj.specification_.is_identity; end
        function value = get.specification_version(obj), value = obj.specification_.specification_version; end
        function value = get.options(obj)
            if numel(obj.stages)>1
                error('MATLAB:Feature:CompositeOptions','A composed Feature has stage-local options; inspect feature.stages.');
            end
            warning('MATLAB:Feature:DeprecatedOptions','Feature.options is a deprecated identity/single-stage convenience; inspect feature.stages.');
            value = struct();
            if ~obj.is_identity, value = obj.stages{1}.options; end
        end
        function varargout = modifier_x0(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_x0(varargin{:});
        end
        function varargout = modifier_affine(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_affine(varargin{:});
        end
        function varargout = modifier_bounds(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_bounds(varargin{:});
        end
        function varargout = modifier_linear_ub(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_linear_ub(varargin{:});
        end
        function varargout = modifier_linear_eq(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_linear_eq(varargin{:});
        end
        function varargout = modifier_fun(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_fun(varargin{:});
        end
        function varargout = modifier_cub(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_cub(varargin{:});
        end
        function varargout = modifier_ceq(obj,varargin)
            kernel = obj.compatibilityKernel();
            [varargout{1:nargout}] = kernel.modifier_ceq(varargin{:});
        end
    end
    methods (Access = private)
        function kernel = compatibilityKernel(obj)
            if numel(obj.stages)>1
                error('MATLAB:Feature:CompositeModifier','Deprecated Feature modifier methods do not execute compositions; use FeaturedProblem.');
            end
            warning('MATLAB:Feature:DeprecatedModifier','Feature modifier methods are deprecated single-stage conveniences; use FeaturedProblem.');
            options = struct();
            if ~obj.is_identity, options = obj.stages{1}.options; end
            kernel = optiprofiler_internal.FeatureKernel(obj.name,options);
        end
    end
    methods (Static)
        function varargout = default_rng(varargin)
            [varargout{1:nargout}] = optiprofiler_internal.FeatureKernel.default_rng(varargin{:});
        end
        function varargout = chebyshevNoiseMap(varargin)
            [varargout{1:nargout}] = optiprofiler_internal.FeatureKernel.chebyshevNoiseMap(varargin{:});
        end
    end
end
