classdef TestProblemReference < matlab.unittest.TestCase
% Feasible reference fact of a problem (the `reference` property of Problem).
%
% The record has exactly four fields: a finite real scalar `merit`, a `kind`
% ('lower_bound', 'optimum', 'best_known' or 'target'; every kind is a claim
% over the FEASIBLE points of the problem), a non-empty `source` and a
% `mapping` token of a closed registry (only 'feasible_objective/1'). It holds
% no point, no function handle and no constraint violation.
%
% What these tests pin:
%
% - the four-field shape, and the rejection of the superseded fields (`fun`,
%   `maxcv`, `point`), of unknown mappings and of non-finite scalars;
% - that a stored record the contract rejects reads as unknown and is never
%   reinterpreted (a file of the superseded layout, a later mapping registry);
% - the propagation table of the features, single and composed, including
%   every subset of the custom option names, uniformly over the four kinds;
% - that a reference changes nothing else: callback counts, histories and
%   whole benchmarks are identical with and without it;
% - what the scalar means: the default merit function returns the objective at
%   every feasible point whatever the initial violation is (the feasible
%   identity), and on a constrained problem the merit of a run may
%   legitimately be below the reference, which is never used as a clamp.
%
% The reference fact must not be confused with two things the code base also
% calls "reference": the reference values of a FeaturedProblem (the scoring
% truth of a trial) and the plain reference of the profiles (a run-history
% minimum). The fact is author/provider metadata and is neither.

    properties (Constant)
        Mapping = 'feasible_objective/1'
        Kinds = {'lower_bound', 'optimum', 'best_known', 'target'}
        CustomKeys = {'mod_affine', 'mod_bounds', 'mod_ceq', 'mod_cub', 'mod_fun', ...
            'mod_linear_eq', 'mod_linear_ub', 'mod_x0'}
        % Written out here, independently of the implementation.
        SafeCustomKeys = {'mod_x0', 'mod_affine'}
    end

    methods (Static)
        function f = shiftedSphere(x)
            f = sum((x(:) - [1; 2]).^2);
        end

        function c = cubFirstCoordinate(x)
            c = x(1) - 3;
        end

        function f = neverEvaluated(~) %#ok<STOUT>
            error('TestProblemReference:evaluated', 'The objective must not be evaluated while loading.');
        end

        function r = record(varargin)
            % record('kind', 'target', ...) overrides fields of a valid record.
            r = struct('merit', 0, 'kind', 'optimum', 'source', 'author', 'mapping', TestProblemReference.Mapping);
            for k = 1:2:numel(varargin)
                r.(varargin{k}) = varargin{k + 1};
            end
        end

        function s = constrainedStruct()
            % Bounds, one linear and one nonlinear constraint, so that every
            % feature has something to act on.
            s = struct('fun', @TestProblemReference.shiftedSphere, 'x0', [0; 0], 'xl', [-5; -5], 'xu', [5; 5], ...
                'aub', [1, 1], 'bub', 10, 'cub', @TestProblemReference.cubFirstCoordinate);
        end

        function problem = constrainedProblem(reference)
            s = TestProblemReference.constrainedStruct();
            if nargin > 0
                s.reference = reference;
            end
            problem = Problem(s);
        end

        function callbacks = customCallbacks()
            % One valid callback per custom option. They are deliberately
            % harmless (the value callbacks return what the predecessor
            % returns, the structure callbacks return the predecessor's
            % structure): the propagation rule reads the option NAMES only,
            % because the framework cannot know what user code does.
            callbacks = struct( ...
                'mod_x0', @(s, p) p.x0 + 0.25, ...
                'mod_affine', @(s, p) deal([2, 0; 1, 1], [0.5; -0.5], inv([2, 0; 1, 1])), ...
                'mod_bounds', @(s, p) deal(p.xl, p.xu), ...
                'mod_linear_ub', @(s, p) deal(p.aub, p.bub), ...
                'mod_linear_eq', @(s, p) deal(p.aeq, p.beq), ...
                'mod_fun', @(x, s, p) p.fun(x), ...
                'mod_cub', @(x, s, p) p.cub(x), ...
                'mod_ceq', @(x, s, p) p.ceq(x));
        end

        function stage = customStage(keys)
            callbacks = TestProblemReference.customCallbacks();
            options = struct();
            for k = 1:numel(keys)
                options.(keys{k}) = callbacks.(keys{k});
            end
            stage = struct('name', 'custom', 'options', options);
        end

        function table = stages()
            % The stage kinds of the propagation table with the expected
            % outcome written as a literal (true: the reference is retained).
            % Non-default options are used on purpose: for a retaining kind
            % no option value matters.
            custom = @TestProblemReference.customStage;
            table = { ...
                'noisy', struct('name', 'noisy', 'options', struct('noise_level', 10)), true; ...
                'noisy_deterministic', struct('name', 'noisy', 'options', struct('noise_mode', 'deterministic')), true; ...
                'truncated', struct('name', 'truncated', 'options', struct('significant_digits', 1)), true; ...
                'random_nan', struct('name', 'random_nan', 'options', struct('nan_rate', 0.9)), true; ...
                'nonquantifiable_constraints', 'nonquantifiable_constraints', true; ...
                'unrelaxable_constraints', struct('name', 'unrelaxable_constraints', 'options', struct( ...
                    'unrelaxable_bounds', true, 'unrelaxable_linear_constraints', true, ...
                    'unrelaxable_nonlinear_constraints', true)), true; ...
                'perturbed_x0', struct('name', 'perturbed_x0', 'options', struct('perturbation_level', 10)), true; ...
                'permuted', 'permuted', true; ...
                'linearly_transformed', struct('name', 'linearly_transformed', 'options', ...
                    struct('rotated', true, 'condition_factor', 5)), true; ...
                'linearly_transformed_unrotated', struct('name', 'linearly_transformed', 'options', ...
                    struct('rotated', false, 'condition_factor', 5)), true; ...
                'quantized_default', 'quantized', false; ...  % ground_truth defaults to true
                'quantized_ground_truth', struct('name', 'quantized', 'options', struct('ground_truth', true)), false; ...
                'quantized_observed', struct('name', 'quantized', 'options', struct('ground_truth', false)), true; ...
                'custom_empty', custom({}), true; ...
                'custom_x0', custom({'mod_x0'}), true; ...
                'custom_affine', custom({'mod_affine'}), true; ...
                'custom_x0_affine', custom({'mod_x0', 'mod_affine'}), true; ...
                'custom_fun', custom({'mod_fun'}), false; ...
                'custom_cub', custom({'mod_cub'}), false; ...
                'custom_ceq', custom({'mod_ceq'}), false; ...
                'custom_bounds', custom({'mod_bounds'}), false; ...
                'custom_linear_ub', custom({'mod_linear_ub'}), false; ...
                'custom_linear_eq', custom({'mod_linear_eq'}), false; ...
                'custom_affine_fun', custom({'mod_affine', 'mod_fun'}), false};
        end

        function [stage, expected] = stageOf(label)
            table = TestProblemReference.stages();
            row = find(strcmp(table(:, 1), label));
            assert(isscalar(row), 'Unknown stage label %s.', label);
            stage = table{row, 2};
            expected = table{row, 3};
        end

        function [M, c] = coordinateMap(featured)
            % Test-local: the map from solver to original coordinates is
            % affine, T(y) = M*y + c, so it is recovered here from n + 1
            % evaluations, without any help from the reference record (which
            % holds no point). The solver coordinates of an original point x
            % are then M \ (x - c).
            n = featured.n;
            c = featured.toOriginalCoordinates(zeros(n, 1));
            M = zeros(n);
            basis = eye(n);
            for j = 1:n
                M(:, j) = featured.toOriginalCoordinates(basis(:, j)) - c;
            end
        end

        function varargout = counted(counter, name, callback, varargin)
            % containers.Map is a handle, so every closure shares one counter.
            counter(name) = counter(name) + 1; %#ok<NASGU>
            [varargout{1:nargout}] = callback(varargin{:});
        end

        function counter = newCounter()
            names = [{'fun', 'cub'}, TestProblemReference.CustomKeys];
            counter = containers.Map(names, num2cell(zeros(1, numel(names))));
        end

        function problem = countedProblem(counter, reference)
            s = TestProblemReference.constrainedStruct();
            s.fun = @(x) TestProblemReference.counted(counter, 'fun', @TestProblemReference.shiftedSphere, x);
            s.cub = @(x) TestProblemReference.counted(counter, 'cub', @TestProblemReference.cubFirstCoordinate, x);
            if ~isempty(reference)
                s.reference = reference;
            end
            problem = Problem(s);
        end

        function stage = countedStage(counter, keys)
            callbacks = TestProblemReference.customCallbacks();
            options = struct();
            for k = 1:numel(keys)
                key = keys{k};
                callback = callbacks.(key);
                if ismember(key, {'mod_fun', 'mod_cub', 'mod_ceq'})
                    options.(key) = @(x, s, p) TestProblemReference.counted(counter, key, callback, x, s, p);
                else
                    options.(key) = @(s, p) TestProblemReference.counted(counter, key, callback, s, p);
                end
            end
            stage = struct('name', 'custom', 'options', options);
        end

        function x = dive(fun, x0, varargin)
            % Visits one slightly infeasible point, then does a short coordinate search.
            fun(x0 - 1.0001);
            x = x0;
            best = fun(x0);
            for i = 1:numel(x0)
                for delta = [0.5, -0.5, 0.25]
                    trial = x;
                    trial(i) = trial(i) + delta;
                    value = fun(trial);
                    if value < best
                        x = trial;
                        best = value;
                    end
                end
            end
        end

        function x = stay(fun, x0, varargin)
            fun(x0);
            x = x0;
        end
    end

    methods (Test)
        % ------------------------------------------------------------ record shape

        function recordHasExactlyFourFields(testCase)
            problem = Problem(struct('fun', @TestProblemReference.shiftedSphere, 'x0', [0; 0], 'reference', ...
                struct('mapping', "feasible_objective/1", 'source', "DOI 10.0/example", 'kind', "best_known", 'merit', int8(3))));
            reference = problem.reference;
            % One canonical field order whatever order the input used, and
            % one canonical type per field.
            testCase.verifyEqual(fieldnames(reference), {'merit'; 'kind'; 'source'; 'mapping'});
            testCase.verifyEqual(reference, struct('merit', 3, 'kind', 'best_known', ...
                'source', 'DOI 10.0/example', 'mapping', 'feasible_objective/1'));
            testCase.verifyClass(reference.merit, 'double');
            testCase.verifyClass(reference.kind, 'char');
            testCase.verifyClass(reference.source, 'char');
            testCase.verifyClass(reference.mapping, 'char');
        end

        function supersededRecordAndTransportApiIsGone(testCase)
            % The superseded layout stored an objective value, a violation and
            % a point, and transported the point through an inverse affine
            % map kept by every stage view. None of this exists any more.
            view = ?optiprofiler_internal.FeatureProblemView;
            testCase.verifyFalse(ismember('affine_inverse', {view.PropertyList.Name}));
            testCase.verifyEmpty(which('optiprofiler_internal.transportProblemReference'));
            testCase.verifyNotEmpty(which('optiprofiler_internal.propagateProblemReference'));
            reference = TestProblemReference.constrainedProblem(TestProblemReference.record()).reference;
            for name = {'fun', 'maxcv', 'point'}
                testCase.verifyFalse(isfield(reference, name{1}), name{1});
            end
        end

        function everyKindIsAcceptedWithTheSameShape(testCase)
            for k = 1:numel(TestProblemReference.Kinds)
                kind = TestProblemReference.Kinds{k};
                reference = TestProblemReference.constrainedProblem(TestProblemReference.record('kind', kind)).reference;
                testCase.verifyEqual(reference.kind, kind);
                testCase.verifyEqual(fieldnames(reference), {'merit'; 'kind'; 'source'; 'mapping'});
            end
        end

        function omittedReferenceIsUnknown(testCase)
            problem = TestProblemReference.constrainedProblem();
            testCase.verifyEmpty(problem.reference);
            testCase.verifyEmpty(TestProblemReference.constrainedProblem([]).reference);
            testCase.verifyEmpty(FeaturedProblem(problem, Feature('noisy'), 5, 0).reference);
            testCase.verifyEmpty(FeaturedProblem(problem, Feature('noisy+permuted'), 5, 0).reference);
        end

        function referenceIsSetAtConstructionOnly(testCase)
            problem = TestProblemReference.constrainedProblem(TestProblemReference.record());
            testCase.verifyError(@() assignReference(problem, []), 'MATLAB:class:SetProhibited');
            testCase.verifyError(@() assignReference(problem, TestProblemReference.record('merit', 1)), ...
                'MATLAB:class:SetProhibited');
            testCase.verifyEqual(problem.reference.merit, 0);
        end

        function supersededFieldsAreRejectedNotReinterpreted(testCase)
            % A record of the superseded layout carries an objective value in
            % `fun`. Reading it as `merit` would be a guess about another
            % contract, so the whole record is refused, even when the four
            % current fields are present as well.
            record = @TestProblemReference.record;  % record(name, value) adds or overrides a field
            legacy = { ...
                struct('fun', 0, 'kind', 'optimum', 'source', 'author'), ...
                struct('fun', 0, 'maxcv', 0, 'kind', 'optimum', 'source', 'author', 'point', [1; 2]), ...
                struct('fun', 6.25, 'maxcv', 0.5, 'kind', 'best_known', 'source', 'test run', 'point', []), ...
                record('fun', 0), record('maxcv', 0), record('point', [1; 2]), record('point', [])};
            for k = 1:numel(legacy)
                testCase.verifyError(@() TestProblemReference.constrainedProblem(legacy{k}), ...
                    'MATLAB:Problem:reference_LegacyField', sprintf('legacy case %d', k));
            end
        end

        function everyFieldIsRequired(testCase)
            % No field has a default. In particular an omitted mapping is not
            % read as 'feasible_objective/1': the author states how the scalar is read.
            for name = {'merit', 'kind', 'source', 'mapping'}
                incomplete = rmfield(TestProblemReference.record(), name{1});
                testCase.verifyError(@() TestProblemReference.constrainedProblem(incomplete), ...
                    'MATLAB:Problem:reference_MissingField', name{1});
            end
        end

        function unknownFieldsAndNakedValuesAreRejected(testCase)
            for name = {'x', 'optimum', 'value', 'merit_fun', 'Merit'}
                extended = TestProblemReference.record(name{1}, 1);
                testCase.verifyError(@() TestProblemReference.constrainedProblem(extended), ...
                    'MATLAB:Problem:reference_UnknownField', name{1});
            end
            % The roadmap forbids a public scalar property: a number alone has
            % no kind, no provenance and no stated reading.
            naked = {1, 0.5, 'optimum', {0, 'optimum', 'author', TestProblemReference.Mapping}, ...
                repmat(TestProblemReference.record(), 1, 2), @sin, true};
            for k = 1:numel(naked)
                testCase.verifyError(@() TestProblemReference.constrainedProblem(naked{k}), ...
                    'MATLAB:Problem:reference_NotStruct', sprintf('naked case %d', k));
            end
        end

        function kindAndSourceAreValidated(testCase)
            cases = { ...
                'kind', 'bogus', 'MATLAB:Problem:reference_kind_Unknown'; ...
                'kind', 'Optimum', 'MATLAB:Problem:reference_kind_Unknown'; ...
                'kind', 'optimum ', 'MATLAB:Problem:reference_kind_Unknown'; ...
                'kind', '', 'MATLAB:Problem:reference_kind_Unknown'; ...
                'kind', 3, 'MATLAB:Problem:reference_kind_NotText'; ...
                'kind', {'optimum'}, 'MATLAB:Problem:reference_kind_NotText'; ...
                'kind', ["optimum", "target"], 'MATLAB:Problem:reference_kind_NotText'; ...
                'kind', string(missing), 'MATLAB:Problem:reference_kind_NotText'; ...
                'source', '', 'MATLAB:Problem:reference_source_Empty'; ...
                'source', '   ', 'MATLAB:Problem:reference_source_Empty'; ...
                'source', "", 'MATLAB:Problem:reference_source_Empty'; ...
                'source', 5, 'MATLAB:Problem:reference_source_NotText'; ...
                'source', ['ab'; 'cd'], 'MATLAB:Problem:reference_source_NotText'};
            for k = 1:size(cases, 1)
                reference = TestProblemReference.record(cases{k, 1}, cases{k, 2});
                testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), cases{k, 3}, ...
                    sprintf('case %d', k));
            end
        end

        % ------------------------------------------------------- mapping registry

        function mappingRegistryIsClosed(testCase)
            % A token is matched exactly. Near misses (another version,
            % another case, surrounding blanks) are different tokens, not
            % spellings of the known one: accepting them would be a guess.
            unknown = {'feasible_objective/2', 'feasible_objective/0', 'feasible_objective', ...
                'feasible_objective/1 ', ' feasible_objective/1', 'Feasible_Objective/1', 'FEASIBLE_OBJECTIVE/1', ...
                'feasible_objective/1.0', 'feasible_objective/01', 'objective/1', 'merit/1', 'default_merit/1', ...
                'identity', ''};
            for k = 1:numel(unknown)
                reference = TestProblemReference.record('mapping', unknown{k});
                testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), ...
                    'MATLAB:Problem:reference_mapping_Unknown', unknown{k});
            end
            % No arbitrary callbacks and no user mappings: a function handle
            % is never a mapping, whatever it computes.
            handles = {@(f, v, v0) f, @sin, @TestProblemReference.shiftedSphere};
            for k = 1:numel(handles)
                reference = TestProblemReference.record('mapping', handles{k});
                testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), ...
                    'MATLAB:Problem:reference_mapping_NotToken', sprintf('handle %d', k));
            end
            other = {1, true, {'feasible_objective/1'}, ["feasible_objective", "1"], ...
                struct('name', 'feasible_objective', 'version', 1)};
            for k = 1:numel(other)
                reference = TestProblemReference.record();
                reference.mapping = other{k};
                testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), ...
                    'MATLAB:Problem:reference_mapping_NotText', sprintf('non-text %d', k));
            end
            % There is no registration function anywhere on the source path.
            source = fileparts(which('Problem'));
            listing = [dir(fullfile(source, '*.m')); dir(fullfile(source, 'private', '*.m')); ...
                dir(fullfile(source, '+optiprofiler_internal', '*.m'))];
            names = lower({listing.name});
            testCase.verifyFalse(any(contains(names, 'mapping')));
        end

        % ------------------------------------------------------- finite validation

        function nonFiniteMeritsAreRejected(testCase)
            values = {NaN, Inf, -Inf, single(NaN), single(Inf), -single(Inf)};
            for k = 1:numel(values)
                for kind = TestProblemReference.Kinds
                    reference = TestProblemReference.record('merit', values{k}, 'kind', kind{1});
                    testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), ...
                        'MATLAB:Problem:reference_merit_NotFinite', sprintf('%d %s', k, kind{1}));
                end
            end
        end

        function nonScalarAndNonRealMeritsAreRejected(testCase)
            % A one-element cell is not a scalar, a logical is not a
            % magnitude, and text is not a number even when it spells one.
            values = {'1', "1", true, false, [], [1, 2], [1; 2], zeros(2), {1}, 1 + 2i, complex(1, 0), ...
                struct('value', 1), @sin};
            for k = 1:numel(values)
                reference = TestProblemReference.record();
                reference.merit = values{k};
                testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), ...
                    'MATLAB:Problem:reference_merit_NotRealScalar', sprintf('case %d', k));
            end
        end

        function integersThatADoubleCannotHoldAreRejected(testCase)
            % Storing 2^53 + 1 as 2^53 would silently change the claim. cast
            % saturates, so intmax is checked explicitly: double(intmax('int64'))
            % is 2^63, which casts back to intmax.
            values = {int64(2)^53 + 1, -int64(2)^53 - 1, intmax('int64'), intmax('uint64'), uint64(2)^53 + 1};
            for k = 1:numel(values)
                reference = TestProblemReference.record('merit', values{k});
                testCase.verifyError(@() TestProblemReference.constrainedProblem(reference), ...
                    'MATLAB:Problem:reference_merit_NotExact', sprintf('case %d', k));
            end
        end

        function finiteRealScalarsAreStoredExactlyAsDouble(testCase)
            values = {0, -3, int64(2)^53, int64(2)^60, intmin('int64'), -0, 1e308, -1e308, 5e-324, ...
                single(0.5), int8(-7), uint16(9), realmax, -realmax};
            stored = [0, -3, 2^53, 2^60, -2^63, 0, 1e308, -1e308, 5e-324, 0.5, -7, 9, realmax, -realmax];
            for k = 1:numel(values)
                reference = TestProblemReference.constrainedProblem( ...
                    TestProblemReference.record('merit', values{k}, 'kind', 'target')).reference;
                testCase.verifyClass(reference.merit, 'double');
                testCase.verifyEqual(reference.merit, stored(k), sprintf('case %d', k));
            end
        end

        function validationNeverEvaluatesTheProblem(testCase)
            counter = TestProblemReference.newCounter();
            without = TestProblemReference.countedProblem(counter, []);
            baseline = [counter('fun'), counter('cub')];
            with_reference = TestProblemReference.countedProblem(counter, TestProblemReference.record());
            % Construction makes exactly the calls it makes without a
            % reference, and none of them is an objective evaluation.
            testCase.verifyEqual([counter('fun'), counter('cub')], 2 * baseline);
            testCase.verifyEqual(counter('fun'), 0);
            testCase.verifyEqual(with_reference.reference.kind, 'optimum');
            testCase.verifyEmpty(without.reference);
            for k = 1:50
                with_reference.reference;
            end
            testCase.verifyEqual([counter('fun'), counter('cub')], 2 * baseline);
            raising = Problem(struct('fun', @TestProblemReference.neverEvaluated, 'x0', [0; 0], ...
                'reference', TestProblemReference.record()));
            % Problem.fun turns the objective's error into a warning and NaN,
            % so the handle is still the raising one: nothing evaluated it.
            testCase.verifyWarning(@() raising.fun([0; 0]), 'TestProblemReference:evaluated');
        end

        function malformedRecordIsRejectedBeforeAnyCallbackIsTouched(testCase)
            % The reference is validated before anything else in the
            % constructor, so a malformed record never causes a callback of
            % the problem to run.
            counter = TestProblemReference.newCounter();
            record = @TestProblemReference.record;  % record(name, value) adds or overrides a field
            malformed = {record('merit', NaN), record('mapping', 'feasible_objective/2'), record('mapping', @sin), ...
                record('kind', 'bogus'), record('source', ''), record('point', [1; 2]), ...
                struct('fun', 0, 'kind', 'optimum', 'source', 'author'), 0};
            for k = 1:numel(malformed)
                try
                    TestProblemReference.countedProblem(counter, malformed{k});
                    testCase.verifyFail(sprintf('malformed case %d was accepted', k));
                catch cause
                    testCase.verifySubstring(cause.identifier, 'MATLAB:Problem:reference_', sprintf('malformed case %d', k));
                end
            end
            testCase.verifyEqual(cell2mat(values(counter)), zeros(1, counter.Count));
        end

        % -------------------------------------------------------- propagation rule

        function ruleMatchesTheTable(testCase)
            % The pure rule; nothing is built or evaluated.
            record = TestProblemReference.record();
            table = TestProblemReference.stages();
            names = {};
            for k = 1:size(table, 1)
                stage = table{k, 2};
                if ischar(stage)
                    name = stage;
                    options = struct();
                else
                    name = stage.name;
                    options = stage.options;
                end
                if strcmp(name, 'quantized') && ~isfield(options, 'ground_truth')
                    options.ground_truth = true;  % the rule sees the defaulted options
                end
                names{end + 1} = name; %#ok<AGROW>
                result = optiprofiler_internal.propagateProblemReference(record, name, options);
                if table{k, 3}
                    testCase.verifyEqual(result, record, table{k, 1});
                else
                    testCase.verifyEmpty(result, table{k, 1});
                end
                % Unknown stays unknown whatever the stage is.
                testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference([], name, options), table{k, 1});
            end
            % A stage kind added later must be classified in the table before
            % it can retain anything: the rule fails closed for other names.
            definitions = optiprofiler_internal.featureDefinitions();
            testCase.verifyEqual(sort(unique(names)), sort(setdiff({definitions.name}, {'plain'})));
            testCase.verifyEqual(optiprofiler_internal.propagateProblemReference(record, 'plain', struct()), record);
            for unknown = {'', 'Noisy', 'rotated', 'mod_affine', 'future_feature'}
                testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference(record, unknown{1}, struct()), ...
                    unknown{1});
            end
        end

        function customWhitelistIsExactlyX0AndAffine(testCase)
            record = TestProblemReference.record();
            definition = optiprofiler_internal.featureDefinitions('custom');
            % The eight custom options, all covered by the subset tests below.
            testCase.verifyEqual(sort(definition.local_keys), TestProblemReference.CustomKeys);
            % A whitelist, not a blacklist: an option name nobody has defined
            % yet makes the reference unknown.
            testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference(record, 'custom', ...
                struct('mod_x0', @sin, 'mod_future', @sin)));
            testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference(record, 'custom', struct('n_runs', 3)));
            testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference(record, 'custom', []));
            % The rule reads names only: callbacks that would raise are never called.
            explode = @(varargin) error('TestProblemReference:called', 'A custom callback was called.');
            testCase.verifyEqual(optiprofiler_internal.propagateProblemReference(record, 'custom', ...
                struct('mod_x0', explode, 'mod_affine', explode)), record);
            testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference(record, 'custom', ...
                struct('mod_fun', explode)));
        end

        function quantizedIsRetainedOnlyForAnExplicitFalse(testCase)
            % Missing or non-logical values count as ground truth (fail closed).
            record = TestProblemReference.record();
            unknown = {struct(), struct('ground_truth', true), struct('ground_truth', 1), struct('ground_truth', 0), ...
                struct('ground_truth', []), struct('ground_truth', 'false'), struct('ground_truth', [false, false]), []};
            for k = 1:numel(unknown)
                testCase.verifyEmpty(optiprofiler_internal.propagateProblemReference(record, 'quantized', unknown{k}), ...
                    sprintf('case %d', k));
            end
            testCase.verifyEqual(optiprofiler_internal.propagateProblemReference(record, 'quantized', ...
                struct('ground_truth', false, 'mesh_size', 0.1)), record);
        end

        % ------------------------------------------------------- propagation table

        function singleStageTableIsUniformOverTheKinds(testCase)
            % The rule is uniform over the kinds: each is a claim over feasible
            % points, so each survives exactly the same stages.
            table = [{'plain', 'plain', true}; TestProblemReference.stages()];
            for i = 1:numel(TestProblemReference.Kinds)
                kind = TestProblemReference.Kinds{i};
                record = TestProblemReference.record('merit', 0.5, 'kind', kind, 'source', ['table:', kind]);
                problem = TestProblemReference.constrainedProblem(record);
                for k = 1:size(table, 1)
                    featured = FeaturedProblem(problem, Feature(table(k, 2)), 10, 3);
                    label = [table{k, 1}, ' ', kind];
                    testCase.verifyNotEqual(featured.execution_strategy, 'composed-views', label);
                    if table{k, 3}
                        % Retained means unchanged: same four fields, nothing transported.
                        testCase.verifyEqual(featured.reference, record, label);
                    else
                        testCase.verifyEmpty(featured.reference, label);
                    end
                end
                testCase.verifyEqual(problem.reference, record);  % the original is untouched
            end
        end

        function everyOrderedPairIsSafeOnlyIfBothStagesAre(testCase)
            table = TestProblemReference.stages();
            record = TestProblemReference.record('kind', 'lower_bound', 'merit', -1);
            problem = TestProblemReference.constrainedProblem(record);
            mismatches = {};
            for i = 1:size(table, 1)
                for j = 1:size(table, 1)
                    featured = FeaturedProblem(problem, Feature({table{i, 2}, table{j, 2}}), 10, 3);
                    assert(strcmp(featured.execution_strategy, 'composed-views'));
                    expected = table{i, 3} && table{j, 3};
                    if expected ~= isequal(featured.reference, record) || expected == isempty(featured.reference)
                        mismatches{end + 1} = [table{i, 1}, '+', table{j, 1}]; %#ok<AGROW>
                    end
                end
            end
            testCase.verifyEmpty(mismatches);
        end

        function longerCompositions(testCase)
            cases = { ...
                {'perturbed_x0', 'noisy', 'permuted', 'linearly_transformed', 'truncated'}, true; ...
                {'custom_x0_affine', 'quantized_observed', 'unrelaxable_constraints', 'random_nan'}, true; ...
                {'quantized_ground_truth', 'noisy', 'permuted'}, false; ...   % unsafe first
                {'noisy', 'custom_bounds', 'permuted'}, false; ...            % unsafe in the middle
                {'noisy', 'permuted', 'custom_fun'}, false; ...               % unsafe last
                % A later safe stage never restores what an earlier stage made unknown.
                {'quantized_default', 'quantized_observed', 'custom_affine'}, false};
            for i = 1:numel(TestProblemReference.Kinds)
                record = TestProblemReference.record('kind', TestProblemReference.Kinds{i});
                problem = TestProblemReference.constrainedProblem(record);
                for k = 1:size(cases, 1)
                    stages = cellfun(@(label) TestProblemReference.stageOf(label), cases{k, 1}, 'UniformOutput', false);
                    featured = FeaturedProblem(problem, Feature(stages), 10, 3);
                    testCase.verifyEqual(featured.execution_strategy, 'composed-views');
                    if cases{k, 2}
                        testCase.verifyEqual(featured.reference, record, strjoin(cases{k, 1}, '+'));
                    else
                        testCase.verifyEmpty(featured.reference, strjoin(cases{k, 1}, '+'));
                    end
                end
            end
        end

        function shorthandOptionsAndPlainTokensReachTheRule(testCase)
            record = TestProblemReference.record();
            problem = TestProblemReference.constrainedProblem(record);
            callbacks = TestProblemReference.customCallbacks();
            testCase.verifyEqual(FeaturedProblem(problem, Feature('plain+noisy+plain'), 10, 3).reference, record);
            testCase.verifyEmpty(FeaturedProblem(problem, Feature('plain+quantized'), 10, 3).reference);
            testCase.verifyEqual(FeaturedProblem(problem, Feature('quantized', struct('ground_truth', false)), 10, 3).reference, record);
            testCase.verifyEqual(FeaturedProblem(problem, Feature('quantized+noisy', struct('ground_truth', false)), 10, 3).reference, record);
            testCase.verifyEqual(FeaturedProblem(problem, Feature('custom', struct('mod_affine', callbacks.mod_affine)), 10, 3).reference, record);
            testCase.verifyEmpty(FeaturedProblem(problem, Feature('custom', struct('mod_fun', callbacks.mod_fun)), 10, 3).reference);
        end

        function allCustomKeySubsets(testCase)
            % All 2^8 subsets of the custom option names, as a single stage and
            % first and last in a composition: retained iff the subset is
            % within {mod_x0, mod_affine}.
            keys = TestProblemReference.CustomKeys;
            record = TestProblemReference.record('kind', 'best_known', 'merit', 0.5);
            problem = TestProblemReference.constrainedProblem(record);
            builders = { ...
                @(stage) Feature({stage}), ...
                @(stage) Feature({stage, 'noisy'}), ...
                @(stage) Feature({'permuted', stage})};
            for b = 1:numel(builders)
                mismatches = {};
                retained = {};
                for mask = 0:(2^numel(keys) - 1)
                    subset = keys(dec2bin(mask, numel(keys)) == '1');
                    expected = all(ismember(subset, TestProblemReference.SafeCustomKeys));
                    stage = TestProblemReference.customStage(subset);
                    reference = FeaturedProblem(problem, builders{b}(stage), 10, 3).reference;
                    if expected ~= isequal(reference, record) || expected == isempty(reference)
                        mismatches{end + 1} = strjoin(subset, ','); %#ok<AGROW>
                    end
                    if ~isempty(reference)
                        retained{end + 1} = strjoin(sort(subset), ','); %#ok<AGROW>
                    end
                    rule = optiprofiler_internal.propagateProblemReference(record, 'custom', stage.options);
                    testCase.assertEqual(~isempty(rule), expected, strjoin(subset, ','));
                end
                testCase.verifyEmpty(mismatches, sprintf('builder %d', b));
                % Exactly the four subsets of {mod_x0, mod_affine} retain the record.
                testCase.verifyEqual(sort(retained), sort({'', 'mod_x0', 'mod_affine', 'mod_affine,mod_x0'}), ...
                    sprintf('builder %d', b));
            end
        end

        % ------------------------------------------------ callback-count invariance

        function aReferenceCostsNoCallbackCall(testCase)
            % The superseded design called the affine modifier once more to
            % transport a point. Nothing is transported now, so building and
            % using a featured problem makes exactly the same calls with a
            % reference as without one, user callbacks included.
            features = { ...
                {'plain'}, {'noisy'}, {'permuted'}, {'linearly_transformed'}, {'perturbed_x0'}, {'quantized'}, ...
                {'unrelaxable_constraints'}, {'noisy', 'permuted'}, {'linearly_transformed', 'quantized', 'noisy'}, ...
                {{'mod_x0', 'mod_affine'}}, {TestProblemReference.CustomKeys}, ...
                {{'mod_affine'}, 'noisy'}, {{'mod_fun', 'mod_cub'}, 'noisy'}};
            for k = 1:numel(features)
                outcome = cell(1, 2);
                counters = cell(1, 2);
                references = {TestProblemReference.record(), []};
                for variant = 1:2
                    counter = TestProblemReference.newCounter();
                    problem = TestProblemReference.countedProblem(counter, references{variant});
                    stages = features{k};
                    for s = 1:numel(stages)
                        if iscell(stages{s})  % a cell of option names stands for a counted custom stage
                            stages{s} = TestProblemReference.countedStage(counter, stages{s});
                        end
                    end
                    featured = FeaturedProblem(problem, Feature(stages), 20, 11);
                    observed = {};
                    for shift = [0, 0.25, -0.5]
                        x = featured.x0 + shift;
                        observed{end + 1} = {featured.fun(x), featured.cub(x), featured.maxcv(x)}; %#ok<AGROW>
                    end
                    before = cell2mat(values(counter));
                    for repeat = 1:100  % reading the fact calls nothing
                        problem.reference;
                        featured.reference;
                    end
                    testCase.verifyEqual(cell2mat(values(counter)), before, sprintf('feature %d', k));
                    counters{variant} = cell2mat(values(counter));
                    outcome{variant} = {featured.x0, observed, featured.fun_hist, featured.maxcv_hist, featured.cub_hist};
                end
                testCase.verifyEqual(counters{1}, counters{2}, sprintf('feature %d', k));
                testCase.verifyGreaterThan(sum(counters{1}), 0);  % the comparison is not vacuous
                testCase.verifyTrue(isequaln(outcome{1}, outcome{2}), sprintf('feature %d', k));
            end
        end

        % ---------------------------------------------------- safe coordinate changes

        function retainedClaimIsTrueForTheFeaturedProblem(testCase)
            % A permutation and a valid affine change of variables keep the
            % reference although the record holds no point. The optimum is
            % known to the test only:
            %
            %   minimize (x1 - 1)^2 + (x2 - 2)^2 + 3  s.t.  x1 + x2 <= 2.5, -5 <= x <= 5.
            %
            % The unconstrained minimizer (1, 2) violates the linear
            % constraint, so the solution is its projection onto the
            % half-plane, x* = (0.75, 1.75), with value 2 * 0.25^2 + 3 = 3.125.
            optimum = 3.125;
            minimizer = [0.75; 1.75];
            record = TestProblemReference.record('merit', optimum, 'source', 'test-local derivation');
            problem = Problem(struct('fun', @(x) TestProblemReference.shiftedSphere(x) + 3, 'x0', [0; 0], ...
                'xl', [-5; -5], 'xu', [5; 5], 'aub', [1, 1], 'bub', 2.5, 'reference', record));
            custom = @TestProblemReference.customStage;
            pipelines = { ...
                {'permuted'}, ...
                {struct('name', 'linearly_transformed', 'options', struct('rotated', true, 'condition_factor', 4))}, ...
                {struct('name', 'linearly_transformed', 'options', struct('rotated', false, 'condition_factor', 4))}, ...
                {custom({'mod_affine'})}, ...
                {custom({'mod_x0', 'mod_affine'})}, ...
                {'permuted', 'noisy'}, ...
                {'linearly_transformed', 'permuted'}, ...
                {custom({'mod_affine'}), 'linearly_transformed', 'truncated'}, ...
                {'perturbed_x0', 'permuted', custom({'mod_affine'}), ...
                    struct('name', 'quantized', 'options', struct('ground_truth', false))}};
            moved = false;
            for k = 1:numel(pipelines)
                for seed = [0, 3, 17]
                    featured = FeaturedProblem(problem, Feature(pipelines{k}), 50, seed);
                    label = sprintf('pipeline %d seed %d', k, seed);
                    testCase.verifyEqual(featured.reference, record, label);

                    % 1. The claimed value is attained at a feasible point of
                    %    the featured problem: the image of x*, found here
                    %    without any stored point.
                    [M, c] = TestProblemReference.coordinateMap(featured);
                    y = M \ (minimizer - c);
                    moved = moved || norm(y - minimizer) > 1e-6;
                    [value, violation] = featured.evaluateTruth(y);
                    testCase.verifyEqual(value, optimum, 'AbsTol', 1e-9, label);
                    testCase.verifyEqual(violation, 0, 'AbsTol', 1e-9, label);

                    % 2. The change of variables is a bijection of the feasible
                    %    set that preserves the truth: at the image of any x,
                    %    the featured truth is the original (objective,
                    %    violation). Hence no feasible point of the featured
                    %    problem is below the reference either.
                    stream = RandStream('mt19937ar', 'Seed', seed);
                    points = -6 + 12 * rand(stream, 2, 300);
                    truth = zeros(2, size(points, 2));     % featured (objective; violation) at the images
                    original = zeros(2, size(points, 2));  % original (objective; violation) at the points
                    for j = 1:size(points, 2)
                        x = points(:, j);
                        [truth(1, j), truth(2, j)] = featured.evaluateTruth(M \ (x - c));
                        original(:, j) = [problem.fun(x); problem.maxcv(x)];
                    end
                    % One qualification per problem: thousands of single
                    % qualifications would dominate the running time.
                    testCase.verifyEqual(truth, original, 'AbsTol', 1e-9, 'RelTol', 1e-9, label);
                    feasible = original(2, :) == 0;
                    testCase.verifyGreaterThanOrEqual(truth(1, feasible), optimum - 1e-9, label);
                    testCase.verifyGreaterThan(nnz(feasible), 50, label);  % the sample really covers the feasible set
                end
            end
            % Guard against a vacuous test: solver coordinates really differ
            % from the original ones, so a stored point could not have been kept.
            testCase.verifyTrue(moved);
        end

        function invalidAffineNeverYieldsAFeaturedProblem(testCase)
            % "Valid affine coordinate change" is enforced where the problem
            % is built: a custom map whose inverse is wrong raises, so it never
            % carries a reference anywhere.
            problem = TestProblemReference.constrainedProblem(TestProblemReference.record());
            broken = struct('name', 'custom', 'options', struct('mod_affine', @(s, p) deal([2, 0; 1, 1], [0; 0], eye(2))));
            testCase.verifyError(@() FeaturedProblem(problem, Feature({broken}), 10, 3), ...
                'MATLAB:Feature:AffineTransformationNotInvertible');
            testCase.verifyError(@() FeaturedProblem(problem, Feature({broken, 'noisy'}), 10, 3), ...
                'MATLAB:Feature:AffineTransformationNotInvertible');
        end

        % ------------------------------------------ default merit: feasible identity

        function defaultMeritIsTheObjectiveAtFeasiblePoints(testCase)
            % 'feasible_objective/1' reads the scalar as an objective value
            % over feasible points. It is comparable with run merits exactly
            % when the merit function maps a feasible point to its objective
            % value, whatever the violation of the run's initial point is.
            current = pwd;
            testCase.addTeardown(@() cd(current));
            cd(fullfile(fileparts(which('Problem')), 'private'));
            initial_violations = [0, 0.2, 50, NaN];
            values = [-1e30, -3.5, -1e-300, -0, 0, 5e-324, 2.5, 1e30, Inf, -Inf];
            for v0 = initial_violations
                for value = values
                    % Exactly the value, not an approximation.
                    testCase.verifyEqual(defaultMerit(value, 0, v0), value, sprintf('f=%g v0=%g', value, v0));
                end
                % The same holds through the array path used for the histories
                % (one initial violation per entry, all equal to v0).
                merits = meritFunCompute(@defaultMerit, values, zeros(size(values)), repmat(v0, size(values)));
                testCase.verifyEqual(merits, values, sprintf('v0=%g', v0));
                for kind = TestProblemReference.Kinds
                    reference = TestProblemReference.record('merit', -7.25, 'kind', kind{1});
                    testCase.verifyEqual(defaultMerit(reference.merit, 0, v0), reference.merit);
                end
            end
            % The identity is a statement about feasible points only. Away
            % from feasibility the default merit tolerates (below v1),
            % penalizes (up to v2) or discards (beyond v2) a point; this is
            % the violation handling that lets run merits fall below the reference.
            tolerated = [1e-10, 1e-10, 5e-9, 1e-10];
            for k = 1:numel(initial_violations)
                v0 = initial_violations(k);
                testCase.verifyEqual(defaultMerit(2, tolerated(k), v0), 2);
                testCase.verifyEqual(defaultMerit(2, 0.05, v0), 2 + 1e5 * (0.05 - tolerated(k)), 'RelTol', 1e-12);
                testCase.verifyEqual(defaultMerit(2, 1e3, v0), Inf);
            end
        end

        function customMeritMustPreserveTheIdentityBeforeAnyConsumerUsesTheReference(testCase)
            % Test-local check that a future consumer has to make before it
            % compares run merits with a record: the probes are the four
            % initial violations above. A merit function without the identity
            % puts run merits in another space than the record.
            current = pwd;
            testCase.addTeardown(@() cd(current));
            cd(fullfile(fileparts(which('Problem')), 'private'));
            shifted = @(f, v, v0) f + 1;
            scaled = @(f, v, v0) f * (1 + v0);
            quadratic = @(f, v, v0) f + 1e3 * v^2;
            testCase.verifyTrue(preservesFeasibleIdentity(@defaultMerit));
            testCase.verifyTrue(preservesFeasibleIdentity(quadratic));  % a custom merit may well have the identity
            testCase.verifyFalse(preservesFeasibleIdentity(shifted));
            testCase.verifyFalse(preservesFeasibleIdentity(scaled));
            % With the shifted merit a feasible optimal run has merit 1, not 0:
            % comparing it with the record would be meaningless, not "close".
            reference = TestProblemReference.record();
            testCase.verifyNotEqual(shifted(reference.merit, 0, 0), reference.merit);
        end

        % ------------------------------------------ constrained merit below reference

        function runMeritsFallBelowTheReferenceAndAreNotClamped(testCase)
            % minimize 1e7 * x subject to -x <= 0: the feasible optimum is 0,
            % at x = 0. An infeasible point has a negative objective, and the
            % default merit does not push it back above 0: below the tolerance
            % v1 it is the objective, and up to v2 the penalty 1e5 * (v - v1)
            % is smaller than the gain 1e7 * v.
            current = pwd;
            testCase.addTeardown(@() cd(current));
            record = TestProblemReference.record('source', 'analytic');
            histories = cell(1, 2);
            for variant = 1:2
                s = struct('fun', @(x) 1e7 * x(1), 'x0', 1, 'cub', @(x) -x(1));
                if variant == 1
                    s.reference = record;
                end
                featured = FeaturedProblem(Problem(s), Feature('plain'), 10, 0);
                % Feasible, optimal, tolerated (v = 5e-11 <= v1 = 1e-10) and
                % penalized (v = 1e-4 in (v1, v2]) points.
                for x = [1, 0, -5e-11, -1e-4]
                    featured.fun(x);
                end
                histories{variant} = {featured.fun_hist, featured.maxcv_hist, featured.fun_init, featured.maxcv_init};
            end
            % The histories are the same with and without a reference.
            testCase.verifyEqual(histories{1}, histories{2});
            testCase.verifyEqual(featured.maxcv_init, 0);
            testCase.verifyEqual(histories{1}{1}, [1e7, 0, -5e-4, -1e3]);
            testCase.verifyEqual(histories{1}{2}, [0, 0, 5e-11, 1e-4]);
            cd(fullfile(fileparts(which('Problem')), 'private'));
            merits = arrayfun(@(f, v) defaultMerit(f, v, 0), histories{1}{1}, histories{1}{2});
            testCase.verifyEqual(merits, [1e7, 0, -5e-4, -1e3 + 1e5 * (1e-4 - 1e-10)], 'RelTol', 1e-12);
            % The feasible points respect the claim; the infeasible ones are
            % below it, in the tolerance zone and in the penalty zone alike.
            testCase.verifyGreaterThanOrEqual(merits(1), record.merit);
            testCase.verifyEqual(merits(2), record.merit);
            testCase.verifyLessThan(merits(3), record.merit);
            testCase.verifyLessThan(merits(4), record.merit);
            % Not even a rigorous lower bound bounds the merits of infeasible
            % points, so no kind is a floor.
            for kind = TestProblemReference.Kinds
                bound = Problem(struct('fun', @(x) 1e7 * x(1), 'x0', 1, 'cub', @(x) -x(1), ...
                    'reference', TestProblemReference.record('kind', kind{1})));
                trial = FeaturedProblem(bound, Feature('plain'), 10, 0);
                trial.fun(-1e-4);
                testCase.verifyLessThan(defaultMerit(trial.fun_hist(end), trial.maxcv_hist(end), trial.maxcv_init), ...
                    trial.reference.merit);
                % Nothing was clamped to the reference, and the run did not
                % move the reference either: it is a fact about the problem,
                % not a minimum of histories.
                testCase.verifyEqual(trial.reference, TestProblemReference.record('kind', kind{1}));
                testCase.verifyEqual(bound.reference, TestProblemReference.record('kind', kind{1}));
            end
        end

        % ------------------------------------------------------ benchmark invariance

        function benchmarkIsIdenticalWithAndWithoutReferences(testCase)
            % Nothing consumes the record: scores, curves and archived
            % histories of a whole benchmark are the same whether the provider
            % states references or not. One provider serves both variants
            % (see fixtures/reference/reference_invariance_load.m), so names, roots and
            % seeds are the same and only the `reference` field differs.
            root = tempname;
            mkdir(root);
            old_cwd = pwd;
            old_path = path;
            old_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            old_switch = getenv('OPTIPROFILER_TEST_STATE_REFERENCES');
            testCase.addTeardown(@() rmdir(root, 's'));
            testCase.addTeardown(@() cd(old_cwd));
            testCase.addTeardown(@() path(old_path));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', old_registry));
            testCase.addTeardown(@() setenv('OPTIPROFILER_TEST_STATE_REFERENCES', old_switch));
            cd(root);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(root, 'registry.mat'));
            fixtures = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'reference');
            registerProblemLibrary(struct('name', 'reference_invariance', 'root', fixtures, ...
                'select_function', 'reference_invariance_select', 'load_function', 'reference_invariance_load'));

            % The provider really serves both variants.
            library = resolveProblemLibrary('reference_invariance');
            setenv('OPTIPROFILER_TEST_STATE_REFERENCES', '1');
            testCase.verifyEqual(library.load('STEEP').reference, TestProblemReference.record('source', 'analytic'));
            testCase.verifyEqual(library.load('SPHERE').reference.kind, 'lower_bound');
            setenv('OPTIPROFILER_TEST_STATE_REFERENCES', '0');
            testCase.verifyEmpty(library.load('STEEP').reference);

            % Scores and curves are compared for every feature. The archive is
            % written only with score_only=false, which also renders the
            % profiles, so the archived histories are compared for one
            % stochastic feature; that keeps the test short.
            features = {'plain', 'noisy', 'permuted+noisy'};
            archived = 'noisy';
            arrays = {'fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs', 'fun_inits', 'maxcv_inits', ...
                'merit_histories', 'merit_outs', 'merit_inits'};
            for k = 1:numel(features)
                outputs = cell(1, 2);
                archives = cell(1, 2);
                for variant = 1:2
                    setenv('OPTIPROFILER_TEST_STATE_REFERENCES', num2str(2 - variant));  % '1' then '0'
                    savepath_variant = fullfile(root, sprintf('f%d-v%d', k, variant));
                    mkdir(savepath_variant);
                    options = struct('plibs', {{'reference_invariance'}}, 'ptype', 'un', 'mindim', 1, 'maxdim', 2, ...
                        'feature_name', features{k}, 'n_runs', 2, 'seed', 5, 'n_jobs', 1, 'max_eval_factor', 8, ...
                        'max_tol_order', 2, 'draw_hist_plots', 'none', 'solver_names', {{'dive', 'stay'}}, ...
                        'benchmark_id', 'invariance', 'savepath', savepath_variant, ...
                        'score_only', ~strcmp(features{k}, archived), 'silent', true);
                    [scores, profile_scores, curves] = benchmark( ...
                        {@TestProblemReference.dive, @TestProblemReference.stay}, options);
                    outputs{variant} = {scores, profile_scores, curves};
                    if strcmp(features{k}, archived)
                        files = dir(fullfile(savepath_variant, '**', 'data_for_loading.mat'));
                        testCase.assertNumElements(files, 1);
                        loaded = load(fullfile(files(1).folder, files(1).name));
                        archives{variant} = loaded.results_plibs{1};
                    end
                end
                testCase.verifyTrue(isequaln(outputs{1}, outputs{2}), features{k});
                testCase.verifyNotEmpty(outputs{1}{3}, features{k});  % the comparison is not vacuous
                if ~strcmp(features{k}, archived)
                    continue
                end
                testCase.verifyEqual(archives{1}.problem_names, archives{2}.problem_names);
                for a = 1:numel(arrays)
                    testCase.verifyTrue(isequaln(archives{1}.(arrays{a}), archives{2}.(arrays{a})), ...
                        [features{k}, ' ', arrays{a}]);
                end
                % The constrained problem was really run below its reference
                % (0), and the archive keeps those merits as they are: never clamped.
                steep = find(strcmp(archives{1}.problem_names, 'STEEP'));
                testCase.assertNumElements(steep, 1);
                merits = archives{1}.merit_histories(steep, :, :, :);
                testCase.verifyLessThan(min(merits(:)), 0, features{k});
            end
        end

        % ---------------------------------------------------- legacy and provider load

        function saveAndLoadRoundTrip(testCase)
            record = TestProblemReference.record('merit', 0.5, 'kind', 'best_known', 'source', 'catalog');
            problem = TestProblemReference.constrainedProblem(record);
            featured = FeaturedProblem(TestProblemReference.constrainedProblem(record), Feature('noisy+permuted'), 10, 3);
            file = [tempname, '.mat'];
            testCase.addTeardown(@() deleteIfPresent(file));
            save(file, 'problem', 'featured');
            lastwarn('');
            loaded = load(file);
            testCase.verifyEmpty(lastwarn());
            testCase.verifyEqual(loaded.problem.reference, record);
            testCase.verifyEqual(loaded.featured.reference, record);
            testCase.verifyEqual(loaded.problem.fun([1; 2]), 0);
        end

        function objectsFromBeforeTheRecordLoadAsUnknown(testCase)
            fixture_root = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'feature-v2');
            lastwarn('');
            loaded = load(fullfile(fixture_root, 'b0-featured-problem.mat'));
            testCase.verifyEmpty(loaded.fp_noisy.reference);
            testCase.verifyEmpty(loaded.fp_noisy.problem.reference);
            testCase.verifyEmpty(loaded.fp_plain.reference);
            testCase.verifyEmpty(loaded.fp_plain.problem.reference);
            testCase.verifyEmpty(lastwarn(), 'Restoring the 1.x fixture must not warn about the new property.');
        end

        function supersededLayoutFileLoadsAsUnknown(testCase)
            % fixtures/reference/legacy-point-record.mat was written by the
            % superseded candidate: `reference` was a stored property holding
            % (fun, maxcv, kind, source, point). Both objects still load and
            % work; their record is refused on every read and its `fun` is
            % never read as a `merit`.
            fixture = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'reference', 'legacy-point-record.mat');
            lastwarn('');
            loaded = load(fixture);
            testCase.verifyEmpty(lastwarn(), 'Loading must not fail or warn; the record is judged when it is read.');
            testCase.verifyClass(loaded.problem, 'Problem');
            testCase.verifyClass(loaded.featured, 'FeaturedProblem');
            testCase.verifyEqual(loaded.problem.name, 'LEGACY');
            testCase.verifyEqual(loaded.problem.fun([1; 2]), 0);
            testCase.verifyEmpty(testCase.verifyWarning(@() loaded.problem.reference, ...
                'MATLAB:Problem:reference_StoredRecordRejected'));
            testCase.verifyEmpty(testCase.verifyWarning(@() loaded.featured.reference, ...
                'MATLAB:Problem:reference_StoredRecordRejected'));
            % A featured problem built from the loaded problem has no reference either.
            rebuilt = testCase.verifyWarning(@() FeaturedProblem(loaded.problem, Feature('permuted'), 5, 0), ...
                'MATLAB:Problem:reference_StoredRecordRejected');
            testCase.verifyEmpty(rebuilt.reference);
        end

        function rejectedStoredRecordsReadAsUnknown(testCase)
            % What a file written under another mapping registry, or a damaged
            % file, leaves in an object (see fixtures/reference/RawReferenceProblem.m).
            % Unknown, never reinterpreted: the problem stays usable, and no
            % field is guessed, defaulted or borrowed from another layout.
            fixtures = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'reference');
            old_path = path;
            testCase.addTeardown(@() path(old_path));
            addpath(fixtures);
            valid = TestProblemReference.record();
            record = @TestProblemReference.record;  % record(name, value) adds or overrides a field
            stored = { ...
                record('mapping', 'feasible_objective/2'), ...   % a token of a later registry
                record('mapping', 'penalized_merit/1'), ...
                record('merit', NaN), record('merit', Inf), ...
                record('kind', 'certified_bound'), record('source', ''), ...
                struct('fun', 0, 'maxcv', 0, 'kind', 'optimum', 'source', 'author', 'point', [1; 2]), ...
                record('point', [1; 2]), rmfield(valid, 'mapping'), 0, 'optimum', {valid}};
            for k = 1:numel(stored)
                problem = RawReferenceProblem(TestProblemReference.constrainedStruct(), stored{k});
                label = sprintf('stored case %d', k);
                testCase.verifyEmpty(testCase.verifyWarning(@() problem.reference, ...
                    'MATLAB:Problem:reference_StoredRecordRejected', label), label);
                testCase.verifyEqual(problem.fun([1; 2]), 0, label);
                file = [tempname, '.mat'];
                save(file, 'problem');
                loaded = load(file);
                delete(file);
                testCase.verifyEmpty(testCase.verifyWarning(@() loaded.problem.reference, ...
                    'MATLAB:Problem:reference_StoredRecordRejected', label), label);
            end
            % A valid stored record is handed out, without a warning.
            problem = RawReferenceProblem(TestProblemReference.constrainedStruct(), valid);
            testCase.verifyEqual(testCase.verifyWarningFree(@() problem.reference), valid);
        end

        function projectX0KeepsTheReference(testCase)
            record = TestProblemReference.record();
            problem = Problem(struct('fun', @TestProblemReference.shiftedSphere, 'x0', [-1; 5], ...
                'xl', [0; 0], 'xu', [3; 3], 'reference', record));
            problem.project_x0();
            testCase.verifyEqual(problem.x0, [0; 3]);
            testCase.verifyEqual(problem.reference, record);
        end

        function providerLoadPreservesTheReferenceWithoutExecutingAnything(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() cleanupProvider(root));
            writeText(fullfile(root, 'reference_toy_select.m'), sprintf([ ...
                'function names = reference_toy_select(~)\n', ...
                '    names = {''REF'', ''PLAIN'', ''NAN'', ''LEGACY'', ''FUTURE'', ''CALLBACK'', ''NAKED''};\n', ...
                'end\n']));
            writeText(fullfile(root, 'reference_toy_never.m'), sprintf([ ...
                'function f = reference_toy_never(~)\n', ...
                '    error(''TestProblemReference:evaluated'', ''The objective was evaluated while loading.'');\n', ...
                'end\n']));
            writeText(fullfile(root, 'reference_toy_load.m'), sprintf([ ...
                'function problem = reference_toy_load(name)\n', ...
                '    s = struct(''fun'', @reference_toy_never, ''x0'', [0; 0], ''name'', name);\n', ...
                '    good = struct(''merit'', 0, ''kind'', ''optimum'', ''source'', ''catalog:reference_toy'', ''mapping'', ''feasible_objective/1'');\n', ...
                '    switch name\n', ...
                '        case ''REF''\n', ...
                '            s.reference = good;\n', ...
                '        case ''NAN''\n', ...
                '            s.reference = good;\n', ...
                '            s.reference.merit = NaN;\n', ...
                '        case ''LEGACY''\n', ...
                '            s.reference = struct(''fun'', 0, ''kind'', ''optimum'', ''source'', ''catalog:reference_toy'', ''point'', [1; 2]);\n', ...
                '        case ''FUTURE''\n', ...
                '            s.reference = good;\n', ...
                '            s.reference.mapping = ''feasible_objective/2'';\n', ...
                '        case ''CALLBACK''\n', ...
                '            s.reference = good;\n', ...
                '            s.reference.mapping = @reference_toy_never;\n', ...
                '        case ''NAKED''\n', ...
                '            s.reference = 0;\n', ...
                '    end\n', ...
                '    problem = Problem(s);\n', ...
                'end\n']));
            registerProblemLibrary(struct('name', 'reference_toy', 'root', root, ...
                'select_function', 'reference_toy_select', 'load_function', 'reference_toy_load'));
            library = resolveProblemLibrary('reference_toy');
            problem = library.load('REF');
            testCase.verifyClass(problem, 'Problem');
            testCase.verifyEqual(problem.reference, struct('merit', 0, 'kind', 'optimum', ...
                'source', 'catalog:reference_toy', 'mapping', 'feasible_objective/1'));  % provenance survives the load
            testCase.verifyEmpty(library.load('PLAIN').reference);
            % A malformed record makes the load fail; the provider is never
            % "helped" by a repaired or reinterpreted record.
            failing = {'NAN', 'MATLAB:Problem:reference_merit_NotFinite'; ...
                'LEGACY', 'MATLAB:Problem:reference_LegacyField'; ...
                'FUTURE', 'MATLAB:Problem:reference_mapping_Unknown'; ...
                'CALLBACK', 'MATLAB:Problem:reference_mapping_NotToken'; ...
                'NAKED', 'MATLAB:Problem:reference_NotStruct'};
            for k = 1:size(failing, 1)
                testCase.verifyError(@() library.load(failing{k, 1}), failing{k, 2}, failing{k, 1});
            end
            % The loaded objective is still the raising callback: nothing evaluated it.
            testCase.verifyWarning(@() problem.fun([0; 0]), 'TestProblemReference:evaluated');
        end
    end
end


function assignReference(problem, value)
    problem.reference = value;
end


function tf = preservesFeasibleIdentity(merit_fun)
    tf = true;
    for value = [-3.5, 0, 2.5]
        for v0 = [0, 0.2, 50, NaN]
            tf = tf && isequal(merit_fun(value, 0, v0), value);
        end
    end
end


function deleteIfPresent(file)
    if exist(file, 'file') == 2
        delete(file);
    end
end


function cleanup = isolateRegistry(testCase)
    original_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    original_pathdef = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF');
    original_startup = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP');
    registry_root = tempname;
    mkdir(registry_root);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(registry_root, 'problem_libraries.mat'));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(registry_root, 'pathdef.m'));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', fullfile(registry_root, 'startup.m'));
    testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', original_registry));
    testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', original_pathdef));
    testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', original_startup));
    testCase.addTeardown(@() removeDirectoryIfPresent(registry_root));
    cleanup = onCleanup(@() rehash);
end


function writeText(filename, contents)
    file_id = fopen(filename, 'w');
    if file_id == -1
        error('Cannot create test provider file %s.', filename);
    end
    cleanup = onCleanup(@() fclose(file_id));
    fprintf(file_id, '%s', contents);
    clear cleanup
end


function removePathIfPresent(root)
    if contains([path, pathsep], [root, pathsep])
        rmpath(root);
    end
end


function removeDirectoryIfPresent(root)
    if exist(root, 'dir') == 7
        rmdir(root, 's');
    end
end


function cleanupProvider(root)
    removePathIfPresent(root);
    removeDirectoryIfPresent(root);
end
