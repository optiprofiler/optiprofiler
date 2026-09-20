classdef Feature < handle
%FEATURE Reusable specification of ordered problem transformations.
%
%   A Feature owns configuration, not an experiment or a live random stream.
%   Each effective stage has its own validated local options. FeaturedProblem
%   builds fresh execution state from the specification for each trial.
%
%   .. rubric:: Construction
%
%   F = Feature(NAME) accepts an atomic char/string name or names joined with
%   '+'. The first name is applied first: 'noisy+truncated' adds noise and then
%   truncates the observed value. Names are trimmed and lowercased.
%
%   F = Feature(NAME, OPTIONS) accepts a scalar struct of local options.
%   F = Feature(NAME, KEY, VALUE, ...) accepts the same options as name/value
%   pairs. A supplied option is broadcast to every declared stage that accepts
%   it; each stage validates it separately. An option accepted by no stage is
%   an error. Repeated stages share these flat supplied values. Option names
%   are case-insensitive; two spellings of one name in the same struct or
%   name/value list are rejected as a duplicate (the configuration would be
%   ambiguous).
%
%   F = Feature(STAGE) accepts a scalar struct with fields 'name' and optional
%   'options' (a scalar struct). F = Feature(STAGES) accepts a nonempty cell
%   array of atomic names or such structs, in order. This form gives repeated
%   stages independent options. Extra flat options, struct arrays, internal
%   stage records, and empty specifications are not accepted. Stage names in
%   structured input must be atomic, not joined with '+'.
%
%   F = Feature(EXISTING_FEATURE) returns the existing canonical Feature
%   without normalizing it again. No extra options are accepted in this form.
%   Feature() without a specification is an error; use Feature('plain').
%
%   'plain' entries are validated before being removed. An all-plain
%   declaration has zero effective stages, while the declaration remains
%   available for inspection. Any number of stages is allowed and names may
%   repeat. For example::
%
%       F = Feature({ ...
%           struct('name','noisy','options',struct('noise_level',1e-2)), ...
%           'plain', ...
%           struct('name','noisy','options',struct('noise_level',1e-4))});
%
%   .. rubric:: Experiment options
%
%   n_runs is never a Feature option, including inside a stage. Supply it to
%   benchmark instead, for example::
%
%       options = struct('feature', F, 'n_runs', 3);
%       scores = benchmark({@solver1, @solver2}, options);
%
%   At the benchmark boundary, 'feature' and 'feature_name' are mutually
%   exclusive. 'feature' accepts a Feature or structured stages, not a bare
%   shorthand string, and cannot be combined with flat local options or load.
%   The experiment layer resolves repetition defaults and the independent
%   plain-reference role; stochastic classification alone does not determine
%   the run count. Experiment-wide controls such as seed are not stage options.
%
%   .. rubric:: Built-in stages and all local options
%
%   The following defaults are local to one stage. A missing local OPTIONS
%   struct uses these defaults. No stage stores n_runs. Numeric options must
%   be finite real scalars and are stored as double; integer-class values must
%   be exactly representable as double. Logical options accept
%   true/false or numeric 0/1 and are stored as logical. Text choices are
%   lowercase and case-sensitive.
%
%   - 'plain': the identity; no local options.
%   - 'perturbed_x0': perturb the initial point. 'distribution' defaults to
%     'spherical'; alternatives are 'gaussian' or a function handle accepting
%     (random_stream, dimension) and returning a perturbation vector.
%     'perturbation_level' is a finite, nonnegative scalar magnitude factor
%     (default 1e-3), scaled by max(1, norm(problem.x0)). MATLAB accepts no
%     vector amplitudes.
%   - 'noisy': modify observed objective and nonlinear-constraint values.
%     'distribution' defaults to 'gaussian'; alternatives are 'uniform' or a
%     function handle accepting (random_stream, output_size).
%     'noise_level' is a finite, nonnegative scalar (default 1e-3).
%     'noise_type' is 'absolute', 'relative', or 'mixed' (default 'mixed').
%     'noise_mode' is 'random' (default) or
%     'deterministic'. In deterministic mode, 'noise_map' is 'chebyshev'
%     (default) or a function handle x -> real scalar. The named noise map is
%     not used in random mode.
%   - 'truncated': truncate observed objective and nonlinear-constraint values.
%     'significant_digits' is the retained number of significant digits, a
%     positive integer (default 6). 'perturbed_trailing_digits' controls
%     randomization of the
%     trailing digits (default false).
%   - 'permuted': randomly permute variables and transport the initial point,
%     bounds and constraints consistently; no local options.
%   - 'linearly_transformed': apply an invertible linear coordinate change.
%     'rotated' controls random rotation (default true). 'condition_factor'
%     is a finite, nonnegative scalar (default 0); the existing
%     transformation has condition number
%     2 ^ sqrt(condition_factor * n / 2) for dimension n >= 2, and 1 for n = 1.
%   - 'random_nan': replace observed objective and nonlinear-constraint values
%     by NaN with probability 'nan_rate', a finite value in [0, 1] (default
%     0.05). These NaN observations are the feature's output; a NaN option
%     value is a configuration error.
%   - 'unrelaxable_constraints': set the observed objective to Inf when a
%     selected category of predecessor constraints is violated.
%     'unrelaxable_bounds' defaults to true;
%     'unrelaxable_linear_constraints' and
%     'unrelaxable_nonlinear_constraints' default to false.
%   - 'nonquantifiable_constraints': return 0 when cub <= 0 or abs(ceq) <= 1e-6,
%     and 1 otherwise, retaining undefined values as NaN; no local options.
%   - 'quantized': evaluate at a point snapped to a mesh. 'mesh_size' is a
%     finite, positive scalar (default 1e-3); 'mesh_type' is 'absolute'
%     (default) or 'relative'. 'ground_truth'
%     defaults to true: the local reference objective and nonlinear
%     constraints also use the snapped point. With false, only observations
%     are snapped. Bounds and linear constraints are assessed at the unsnapped
%     point in that stage's coordinate system.
%   - 'custom': user-supplied modifier function handles, listed below. Omitted
%     modifiers use the runtime's default behavior; a supplied mod_affine still
%     changes coordinates and transports omitted structural components.
%
%   Custom local options have these signatures; random_stream is supplied by
%   the runtime, and problem is the stage's immediate predecessor Problem::
%
%       mod_x0:        (random_stream, problem) -> modified_x0
%       mod_affine:    (random_stream, problem) -> (A, b, inverse_A)
%       mod_bounds:    (random_stream, problem) -> (modified_xl, modified_xu)
%       mod_linear_ub: (random_stream, problem) -> (modified_aub, modified_bub)
%       mod_linear_eq: (random_stream, problem) -> (modified_aeq, modified_beq)
%       mod_fun:       (x, random_stream, problem) -> modified_fun
%       mod_cub:       (x, random_stream, problem) -> modified_cub
%       mod_ceq:       (x, random_stream, problem) -> modified_ceq
%
%   mod_affine uses the coordinate map A*x+b and its inverse. It is asked once
%   per problem and seed, and the initial point, the bounds, the linear
%   constraints and every evaluation use that one answer. The triple is
%   validated when the problem is built: real, finite arrays of matching sizes,
%   norm(abs(inverse_A)*abs(A),inf) < 1/eps, and inverse_A inverts A from both
%   sides, norm(max(abs(A*inverse_A-eye(n)) -
%   64*n*eps*abs(A)*abs(inverse_A),0),'fro') <= 1e-8*n and the same for
%   inverse_A*A (each entry is forgiven the rounding of the products it is
%   summed from, and nothing more); otherwise an error is raised. Integer and
%   single precision arrays are converted to double precision, which rounds to
%   nearest; the rounded triple is the one that is validated and used. If A is
%   exactly diagonal and diag(inverse_A) is its reciprocal to roundoff, the
%   bounds stay bounds; otherwise every finite bound is posed as a linear
%   constraint, so an off-diagonal entry of A is never ignored, however small.
%   Nothing finite is lost: a bound, coefficient or right-hand side that
%   overflows, or that underflows (a nonzero bound below realmin, or an entry
%   of a transported row or right-hand side whose terms are all below it),
%   raises an error. The initial point is inverse_A*(x0-b) only if A maps it
%   back to x0 to the rounding of that evaluation, in every component;
%   otherwise A*y = x0-b is solved, and an error is raised if that point is not
%   mapped back either. Since mod_linear_ub and mod_linear_eq replace the
%   linear constraints, supplying one of them with an A that is not diagonal
%   raises an error if the problem has such bounds, unless mod_bounds is
%   supplied as well. mod_bounds replaces the bounds and nothing else: under an
%   A for which the bounds stay bounds the bounds of the problem are replaced,
%   under any other A they are not re-added as generated rows; explicit linear
%   constraints are still transported. The derivative methods of a featured
%   problem follow the same transformation by the chain rule (see
%   FeaturedProblem). Within an
%   observation callback, querying problem.fun/cub/ceq serves that predecessor
%   again; these calls are not silently treated as reference reads.
%
%   .. rubric:: Inspection, transport and compatibility
%
%   Read-only properties are 'stages', 'declared', 'name', 'declared_name',
%   'is_stochastic', 'is_identity' and 'specification_version'. 'stages' is a
%   cell array of effective records with name, occurrence (same-kind zero-based
%   index), identity (for example 'noisy#0'), literal code and normalized local
%   options. The effective 'name' joins those stages, or is 'plain' for identity.
%   'is_identity' tests for zero effective stages. 'is_stochastic' reports the
%   stages' established classification (custom is classified as stochastic),
%   not the number of actual executions or a guarantee of independent draws.
%
%   'declared' retains a route ('feature_name' or 'feature') and entries with
%   supplied name/options before defaults, including plain entries;
%   'declared_name' joins those declared names. An imported historical Feature
%   with no recorded declaration has declared=[] and declared_name=''. The
%   specification's declaration route is distinct from the benchmark keyword
%   through which an existing Feature is later supplied.
%
%   Returned records and option structs are value copies. Native function
%   handles and handle-valued user state remain references; Feature does not
%   promise deep immutability of a callback's captured state. Normalization
%   does not invoke callbacks; their outputs are checked during execution.
%
%   Current native save/load transports versioned, resolved effective options
%   and the separately retained declaration, not runtime counters or streams.
%   The known old stored name/options layout loads into a LegacyFeatureEnvelope.
%   optiprofiler_internal.importLegacyFeature converts that envelope to a clean
%   Feature plus a separate retained experiment request; an old resolved count
%   is not proof that the old caller requested it explicitly. For saved options,
%   loadBenchmarkOptions prepares a fresh benchmark rather than replotting or
%   resuming a live runtime. Native MAT loading is trusted input and can execute
%   code; original callback/class dependencies must be available. JSON callback
%   descriptions are not executable reconstruction recipes.
%
%   'options' and modifier_x0/affine/bounds/linear_ub/linear_eq/fun/cub/ceq are
%   deprecated identity/single-effective-stage conveniences. They warn, and
%   accessing them for multiple effective stages is an error. Each modifier
%   call creates a fresh numerical kernel; the engine does not use these
%   conveniences to execute a pipeline. Use FeaturedProblem for execution. In
%   particular, a kernel keeps the transformation of mod_affine for one problem
%   and seed, so that the initial point, the bounds, the linear constraints,
%   the derivatives and every evaluation of a FeaturedProblem share one map,
%   whereas each of these conveniences asks mod_affine again: a callback with
%   a state can give two of them two different maps.
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
