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
        function saved = saveobj(obj)
            % Trusted native transport preserves callbacks, not live runtime
            % state, run counts, or diagnostic JSON callback descriptions.
            saved = struct('schema','matlab-feature-native-v2', ...
                'specification',obj.specification_);
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
        function obj = loadobj(saved)
            if strcmp(class(saved),'Feature') && isscalar(saved)
                if ~strcmp(saved.specification_version,'matlab-feature-spec-v2')
                    invalidNative('The canonical specification version is unsupported.');
                end
                obj = saved;
                return;
            end
            if hasExactFields(saved,{'name','options'})
                obj = optiprofiler_internal.LegacyFeatureEnvelope(saved);
                return;
            end
            if ~hasExactFields(saved,{'schema','specification'}) || ...
                    ~isNativeText(saved.schema) || ~strcmp(saved.schema,'matlab-feature-native-v2')
                invalidNative('The native Feature schema is unsupported.');
            end
            state = saved.specification;
            required = {'stages','declared','declared_name','specification_version'};
            allowed = [required,{'name','is_stochastic','is_identity'}];
            if ~isa(state,'struct') || ~isscalar(state) || ...
                    ~all(isfield(state,required)) || ~all(ismember(fieldnames(state),allowed)) || ...
                    ~isNativeText(state.specification_version) || ...
                    ~strcmp(state.specification_version,'matlab-feature-spec-v2') || ...
                    ~isa(state.stages,'cell') || ~(isrow(state.stages) || isempty(state.stages))
                invalidNative('The native canonical specification layout is unsupported.');
            end
            entries = cell(1,numel(state.stages));
            for k = 1:numel(state.stages)
                stage = state.stages{k};
                if ~hasExactFields(stage,{'name','occurrence','identity','code','options'}) || ...
                        ~isNativeText(stage.name) || ~isNativeText(stage.identity) || ...
                        ~isa(stage.options,'struct') || ~isscalar(stage.options) || ...
                        ~nativeNumber(stage.occurrence) || ~nativeNumber(stage.code)
                    invalidNative('A native stage record is malformed.');
                end
                entries{k} = struct('name',stage.name,'options',stage.options);
            end
            % Exactly one construction normalizes local options. The saved
            % derived identity fields are checked, never trusted as authority.
            if isempty(entries), obj = Feature('plain'); else, obj = Feature(entries); end
            canonical = obj.stages;
            if numel(canonical) ~= numel(state.stages)
                invalidNative('Native effective stages must not contain plain.');
            end
            for k = 1:numel(canonical)
                expected = canonical{k}; retained = state.stages{k};
                if ~strcmp(expected.name,retained.name) || ~strcmp(expected.identity,retained.identity) || ...
                        expected.occurrence ~= retained.occurrence || expected.code ~= retained.code
                    invalidNative('Native stage identity, occurrence, or code is inconsistent.');
                end
            end
            if isfield(state,'name') && (~isNativeText(state.name) || ~strcmp(state.name,obj.name))
                invalidNative('The native effective feature name is inconsistent.');
            end
            for field = {'is_stochastic','is_identity'}
                key = field{1};
                if isfield(state,key) && (~isa(state.(key),'logical') || ~isscalar(state.(key)) || state.(key) ~= obj.(key))
                    invalidNative('A native derived feature flag is inconsistent.');
                end
            end
            validateNativeDeclaration(state.declared,state.declared_name,canonical);
            obj.specification_.declared = state.declared;
            obj.specification_.declared_name = state.declared_name;
        end
        function varargout = default_rng(varargin)
            [varargout{1:nargout}] = optiprofiler_internal.FeatureKernel.default_rng(varargin{:});
        end
        function varargout = chebyshevNoiseMap(varargin)
            [varargout{1:nargout}] = optiprofiler_internal.FeatureKernel.chebyshevNoiseMap(varargin{:});
        end
    end
end

function valid = hasExactFields(value,expected)
    valid = isa(value,'struct') && isscalar(value) && isequal(sort(fieldnames(value)),sort(expected(:)));
end
function valid = isNativeText(value)
    valid = isa(value,'char') && (isrow(value) || isempty(value));
end
function valid = nativeNumber(value)
    valid = builtin('isnumeric',value) && isreal(value) && isscalar(value) && isfinite(value);
end
function validateNativeDeclaration(declared,declared_name,stages)
    if ~isNativeText(declared_name), invalidNative('The native declared name must be text.'); end
    if isa(declared,'double') && isempty(declared) && isequal(size(declared),[0,0])
        if ~isempty(declared_name), invalidNative('An unknown declaration cannot have a declared name.'); end
        return;
    end
    if ~hasExactFields(declared,{'route','entries'}) || ~isNativeText(declared.route) || ...
            ~ismember(declared.route,{'feature_name','feature'}) || ~isa(declared.entries,'cell') || ...
            ~isrow(declared.entries) || isempty(declared.entries)
        invalidNative('The native declaration layout is malformed.');
    end
    names = cell(1,numel(declared.entries));
    for k = 1:numel(declared.entries)
        entry = declared.entries{k};
        if ~hasExactFields(entry,{'name','options'}) || ~isNativeText(entry.name) || ...
                ~isa(entry.options,'struct') || ~isscalar(entry.options)
            invalidNative('A native declaration entry is malformed.');
        end
        definition = optiprofiler_internal.featureDefinitions(entry.name);
        if ~strcmp(definition.name,entry.name) || ~all(ismember(fieldnames(entry.options),definition.local_keys))
            invalidNative('A native declaration contains nonlocal options.');
        end
        names{k} = entry.name;
    end
    % Filtering one plain entry yields 0-by-0; compare canonical row sequences.
    effective = reshape(names(~strcmp(names,'plain')),1,[]);
    expected = reshape(cellfun(@(s) s.name,stages,'UniformOutput',false),1,[]);
    if ~strcmp(strjoin(names,'+'),declared_name) || ~isequal(effective,expected)
        invalidNative('The retained declaration is inconsistent with the effective stage order.');
    end
end
function invalidNative(message)
    error('MATLAB:Feature:InvalidNativeState','%s',message);
end
