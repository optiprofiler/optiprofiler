classdef TestFeatureValueValidation < matlab.unittest.TestCase
% Stage values are validated and canonicalized before any runtime state exists.
% e7d5341 accepted the rejected rows below: any perturbation_level, NaN or Inf
% magnitudes, and integer or single classes that later corrupted or crashed
% the kernels. It also stored numeric 0/1 logical options as double, which
% changed their effective-stage provenance.

    properties (TestParameter)
        option = struct( ...
            'perturbation_level', {{'perturbed_x0', 'perturbation_level', 'MATLAB:Feature:perturbation_level_InvalidInput'}}, ...
            'noise_level', {{'noisy', 'noise_level', 'MATLAB:Feature:noise_level_NotPositive'}}, ...
            'condition_factor', {{'linearly_transformed', 'condition_factor', 'MATLAB:Feature:condition_factor_InvalidInput'}}, ...
            'mesh_size', {{'quantized', 'mesh_size', 'MATLAB:Feature:mesh_size_NotPositive'}}, ...
            'nan_rate', {{'random_nan', 'nan_rate', 'MATLAB:Feature:nan_rate_NotBetween_0_1'}});
        invalid = struct('nan_value', NaN, 'inf_value', Inf, 'minus_inf', -Inf, 'negative', -0.1, ...
            'char_text', 'bad', 'string_text', "bad", 'empty_value', [], 'row_vector', [0.1 0.2], ...
            'column_vector', [0.1; 0.2], 'complex_value', 0.1i, 'logical_value', true, 'cell_value', {{0.1}});
        route = {'shorthand', 'structured'};
        logicalOption = struct('rotated', {{'linearly_transformed', 'rotated'}}, ...
            'perturbed_trailing_digits', {{'truncated', 'perturbed_trailing_digits'}}, ...
            'ground_truth', {{'quantized', 'ground_truth'}}, ...
            'unrelaxable_bounds', {{'unrelaxable_constraints', 'unrelaxable_bounds'}});
    end

    methods (Test)
        function invalidMagnitudesAreRejected(testCase, option, invalid, route)
            testCase.verifyError(@() TestFeatureValueValidation.make(option{1}, option{2}, invalid, route), option{3});
        end

        function zeroAndBoundPolicy(testCase, option)
            [name, key, identifier] = option{:};
            if strcmp(key, 'mesh_size')
                testCase.verifyError(@() Feature(name, key, 0), identifier);
            else
                feature = Feature(name, key, 0);
                testCase.verifyEqual(feature.stages{1}.options.(key), 0);
            end
            if strcmp(key, 'nan_rate')
                testCase.verifyError(@() Feature(name, key, 1.5), identifier);
                feature = Feature(name, key, 1);
                testCase.verifyEqual(feature.stages{1}.options.(key), 1);
            end
        end

        function significantDigitsIsAPositiveInteger(testCase)
            for value = {NaN, Inf, 0, -1, 2.5, true, 'six', [], [2 3]}
                testCase.verifyError(@() Feature('truncated', 'significant_digits', value{1}), ...
                    'MATLAB:Feature:significant_digits_NotPositiveInteger');
            end
            feature = Feature('truncated', 'significant_digits', int32(3));
            testCase.verifyClass(feature.stages{1}.options.significant_digits, 'double');
            testCase.verifyEqual(feature.stages{1}.options.significant_digits, 3);
        end

        function numericClassesAreStoredAsDouble(testCase)
            cases = {'perturbed_x0', 'perturbation_level', uint8(1); 'noisy', 'noise_level', int32(1); ...
                'noisy', 'noise_level', single(0.25); 'linearly_transformed', 'condition_factor', int32(2); ...
                'quantized', 'mesh_size', single(0.5); 'random_nan', 'nan_rate', single(0.25)};
            for i = 1:size(cases, 1)
                feature = Feature(cases{i, 1}, cases{i, 2}, cases{i, 3});
                stored = feature.stages{1}.options.(cases{i, 2});
                testCase.verifyClass(stored, 'double');
                testCase.verifyEqual(stored, double(cases{i, 3}));
            end
        end

        function kernelsComputeInDouble(testCase)
            constant = Problem(struct('fun', @(x) 0.3, 'x0', 0));
            integer = FeaturedProblem(constant, Feature('noisy', struct('noise_level', int32(1), 'noise_type', 'absolute')), 10, 17);
            reference = FeaturedProblem(constant, Feature('noisy', struct('noise_level', 1, 'noise_type', 'absolute')), 10, 17);
            value = integer.fun(0);
            testCase.verifyClass(value, 'double');
            testCase.verifyEqual(value, reference.fun(0));
            testCase.verifyNotEqual(value, 0.3);

            witness = Problem(struct('fun', @(x) double(isa(x, 'double')), 'x0', 0));
            quantized = FeaturedProblem(witness, Feature('quantized', struct('mesh_size', int32(1))), 10, 17);
            testCase.verifyEqual(quantized.fun(2.6), 1, 'The user callback must receive a double point.');

            twoD = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2]));
            a = FeaturedProblem(twoD, Feature('perturbed_x0', 'perturbation_level', int32(1)), 10, 17);
            b = FeaturedProblem(twoD, Feature('perturbed_x0', 'perturbation_level', 1), 10, 17);
            testCase.verifyClass(a.x0, 'double');
            testCase.verifyEqual(a.x0, b.x0);
            a = FeaturedProblem(twoD, Feature('linearly_transformed', struct('condition_factor', int32(2), 'rotated', false)), 10, 17);
            b = FeaturedProblem(twoD, Feature('linearly_transformed', struct('condition_factor', 2, 'rotated', false)), 10, 17);
            testCase.verifyEqual(a.x0, b.x0);

            trunc = Problem(struct('fun', @(x) 1.23456789, 'x0', 0));
            options = struct('perturbed_trailing_digits', true, 'significant_digits', int32(3));
            a = FeaturedProblem(trunc, Feature('truncated', options), 10, 17);
            options.significant_digits = 3;
            b = FeaturedProblem(trunc, Feature('truncated', options), 10, 17);
            value = a.fun(0);
            testCase.verifyEqual(value, b.fun(0));
            testCase.verifyNotEqual(value, round(1.23456789, 2), 'The trailing digits must stay perturbed.');
        end

        function numericLogicalsAreStoredAsLogical(testCase, logicalOption)
            [name, key] = logicalOption{:};
            numeric = Feature(name, key, 1);
            logical_ = Feature(name, key, true);
            testCase.verifyClass(numeric.stages{1}.options.(key), 'logical');
            testCase.verifyEqual(numeric.stages{1}.options.(key), true);
            % Effective stages are canonical; the declaration keeps the raw input.
            a = optiprofiler_internal.featureProvenance(numeric, [], struct());
            b = optiprofiler_internal.featureProvenance(logical_, [], struct());
            testCase.verifyEqual(jsonencode(a.feature.stages), jsonencode(b.feature.stages));
        end

        function perturbationLevelIsScalarOnly(testCase)
            % A row amplitude added one inner-product value to every coordinate
            % and a column amplitude failed at run time. Vectors are Python-only.
            for value = {[0.1 0.2 0.3], [0.1; 0.2; 0.3], zeros(0, 1)}
                testCase.verifyError(@() Feature('perturbed_x0', 'perturbation_level', value{1}), ...
                    'MATLAB:Feature:perturbation_level_InvalidInput');
            end
        end

        function randomNanObservationsAreNotConfiguration(testCase)
            fp = FeaturedProblem(Problem(struct('fun', @(x) x(1), 'x0', 0)), Feature('random_nan', 'nan_rate', 1), 20, 17);
            values = arrayfun(@(x) fp.fun(x), linspace(-1, 1, 10));
            testCase.verifyTrue(all(isnan(values)));
        end

        function meshTypeStaysCaseSensitive(testCase)
            testCase.verifyError(@() Feature('quantized', 'mesh_type', 'RELATIVE'), 'MATLAB:Feature:mesh_type_InvalidInput');
            p = Problem(struct('x0', 2.6, 'fun', @(x) x(1)));
            relative = FeaturedProblem(p, Feature('quantized', struct('mesh_type', 'relative', 'mesh_size', 0.5, 'ground_truth', false)), 10, 17);
            absolute = FeaturedProblem(p, Feature('quantized', struct('mesh_type', 'absolute', 'mesh_size', 0.5, 'ground_truth', false)), 10, 17);
            testCase.verifyEqual(relative.fun(2.6), 2.6, 'AbsTol', 1e-15);
            testCase.verifyEqual(absolute.fun(2.6), 2.5, 'AbsTol', 1e-15);
        end

        function invalidHistoricalStateFailsExplicitly(testCase)
            stage = struct('name', 'perturbed_x0', 'occurrence', 0, 'identity', 'perturbed_x0#0', 'code', 1, ...
                'options', struct('distribution', 'spherical', 'perturbation_level', NaN));
            testCase.verifyError(@() Feature.loadobj(TestFeatureValueValidation.nativeState(stage)), ...
                'MATLAB:Feature:perturbation_level_InvalidInput');
            envelope = optiprofiler_internal.LegacyFeatureEnvelope(struct('name', 'noisy', 'options', struct('noise_level', Inf)));
            testCase.verifyError(@() optiprofiler_internal.importLegacyFeature(envelope), 'MATLAB:Feature:noise_level_NotPositive');
        end

        function historicalNumericLogicalIsCanonicalizedOnLoad(testCase)
            stage = struct('name', 'linearly_transformed', 'occurrence', 0, 'identity', 'linearly_transformed#0', 'code', 5, ...
                'options', struct('rotated', 1, 'condition_factor', 0));
            feature = Feature.loadobj(TestFeatureValueValidation.nativeState(stage));
            testCase.verifyClass(feature.stages{1}.options.rotated, 'logical');
            testCase.verifyTrue(feature.stages{1}.options.rotated);
        end

        function integerValuesMustBeExactDoubles(testCase)
            % Storing as double must not round an integer-class input.
            testCase.verifyError(@() Feature('noisy', 'noise_level', int64(2)^53 + 1), ...
                'MATLAB:Feature:noise_level_NotPositive');
            testCase.verifyError(@() Feature('truncated', 'significant_digits', uint64(2)^60 + 1), ...
                'MATLAB:Feature:significant_digits_NotPositiveInteger');
            feature = Feature('noisy', 'noise_level', int64(2)^53);
            testCase.verifyEqual(feature.stages{1}.options.noise_level, 2^53);
            feature = Feature('noisy', 'noise_level', single(0.1));
            testCase.verifyEqual(feature.stages{1}.options.noise_level, double(single(0.1)));
        end

        function saturatingIntegerCastsAreNotExact(testCase)
            % intmax rounds up when converted to double, then saturates back
            % on the reverse cast. That round trip is not proof of exactness.
            options = {'noisy', 'noise_level', 'noise_level_NotPositive'; ...
                'perturbed_x0', 'perturbation_level', 'perturbation_level_InvalidInput'; ...
                'quantized', 'mesh_size', 'mesh_size_NotPositive'; ...
                'linearly_transformed', 'condition_factor', 'condition_factor_InvalidInput'; ...
                'truncated', 'significant_digits', 'significant_digits_NotPositiveInteger'};
            for i = 1:size(options, 1)
                for value = {intmax('int64'), intmax('uint64')}
                    for input_route = {'shorthand', 'structured'}
                        testCase.verifyError(@() TestFeatureValueValidation.make( ...
                            options{i, 1}, options{i, 2}, value{1}, input_route{1}), ...
                            ['MATLAB:Feature:', options{i, 3}]);
                    end
                end
            end
            % Exact integers above flintmax remain valid; do not impose a
            % blanket flintmax ceiling to compensate for the saturating cast.
            for value = {int64(2)^53, uint64(2)^63, ...
                    intmax('int64')-int64(1023), intmax('uint64')-uint64(2047), intmax('uint32')}
                feature = Feature('noisy', 'noise_level', value{1});
                testCase.verifyEqual(feature.stages{1}.options.noise_level, double(value{1}));
            end
        end
    end

    methods (Static)
        function f = make(name, key, value, route)
            if strcmp(route, 'shorthand')
                f = Feature(name, key, value);
            else
                options = struct();
                options.(key) = value;
                f = Feature(struct('name', name, 'options', options));
            end
        end

        function state = nativeState(stage)
            specification = struct('stages', {{stage}}, 'declared', [], 'declared_name', '', ...
                'specification_version', 'matlab-feature-spec-v2');
            state = struct('schema', 'matlab-feature-native-v2', 'specification', specification);
        end
    end
end
