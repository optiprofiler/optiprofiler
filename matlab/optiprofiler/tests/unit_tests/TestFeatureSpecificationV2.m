classdef TestFeatureSpecificationV2 < matlab.unittest.TestCase
% Public construction and inspection of a reusable local-only specification.
    methods (Test)
        function repeatedStagesAreLocalAndReusable(testCase)
            supplied = {struct('name','noisy','options',struct('noise_level',.01)), ...
                'plain', struct('name','noisy','options',struct('noise_level',.0001))};
            feature = Feature(supplied);
            reused = Feature(feature);
            testCase.verifyEqual(reused.name, 'noisy+noisy');
            stages = reused.stages;
            testCase.verifyEqual(numel(stages), 2);
            testCase.verifyEqual(stages{1}.identity, 'noisy#0');
            testCase.verifyEqual(stages{2}.identity, 'noisy#1');
            testCase.verifyEqual(stages{1}.options.noise_level, .01);
            testCase.verifyEqual(stages{2}.options.noise_level, .0001);
            testCase.verifyFalse(isfield(stages{1}.options, 'n_runs'));
            testCase.verifyFalse(isfield(stages{2}.options, 'n_runs'));
            testCase.verifyEqual(reused.declared_name, 'noisy+plain+noisy');
        end
        function allKindsAndShorthandHaveCanonicalLocalRecords(testCase)
            names = {'perturbed_x0','noisy','truncated','permuted','linearly_transformed', ...
                'random_nan','unrelaxable_constraints','nonquantifiable_constraints','quantized','custom'};
            for k = 1:numel(names)
                feature = Feature(names{k});
                stage = feature.stages{1};
                testCase.verifyEqual(stage.name, names{k});
                testCase.verifyEqual(stage.code, k);
                testCase.verifyFalse(isfield(stage.options,'n_runs'));
            end
            shorthand = Feature(' plain + noisy + perturbed_x0 + noisy ', ...
                struct('distribution','gaussian','noise_level',.02));
            explicit = Feature({'plain', ...
                struct('name','noisy','options',struct('distribution','gaussian','noise_level',.02)), ...
                struct('name','perturbed_x0','options',struct('distribution','gaussian')), ...
                struct('name','noisy','options',struct('distribution','gaussian','noise_level',.02))});
            testCase.verifyEqual(shorthand.stages, explicit.stages);
            testCase.verifyEqual(shorthand.declared.route, 'feature_name');
            testCase.verifyEqual(explicit.declared.route, 'feature');
            testCase.verifyEqual(shorthand.declared_name, 'plain+noisy+perturbed_x0+noisy');
            identity = Feature('plain+plain');
            testCase.verifyEmpty(identity.stages);
            testCase.verifyTrue(identity.is_identity);
            testCase.verifyFalse(identity.is_stochastic);
        end
        function ownedInspectionAndDeprecatedNumericsRemainUsable(testCase)
            feature = Feature('plain');
            testCase.verifyWarning(@() feature.options, 'MATLAB:Feature:DeprecatedOptions');
            state = warning; restore = onCleanup(@() warning(state));
            warning('off','MATLAB:Feature:DeprecatedOptions');
            warning('off','MATLAB:Feature:DeprecatedModifier');
            testCase.verifyEqual(feature.options, struct());
            problem = Problem(struct('fun',@(x) x'*x,'x0',[1;2]));
            testCase.verifyEqual(feature.modifier_fun([1;2],17,problem,0),5);
            stream = Feature.default_rng(17);
            testCase.verifyClass(stream,'RandStream');
            testCase.verifyError(@() Feature.default_rng('a'),'MATLAB:Feature:SeedNotEvenReal');
            supplied = struct('noise_level',.01);
            feature = Feature('noisy',supplied);
            supplied.noise_level = 5;
            observed = feature.stages; observed{1}.options.noise_level = 6;
            declared = feature.declared; declared.entries{1}.options.noise_level = 7;
            testCase.verifyEqual(feature.stages{1}.options.noise_level,.01);
            testCase.verifyEqual(feature.declared.entries{1}.options.noise_level,.01);
            testCase.verifyError(@() Feature('plain','n_runs',3),'MATLAB:Feature:ExperimentOption');
            testCase.verifyError(@() Feature(struct('name','plain','options',struct('noise_level',.1))), ...
                'MATLAB:Feature:InvalidOptionForFeature');
            testCase.verifyError(@() Feature({'noisy'},'noise_level',.1),'MATLAB:Feature:StructuredOverrides');
            multiple = Feature('noisy+noisy');
            testCase.verifyError(@() multiple.options,'MATLAB:Feature:CompositeOptions');
            testCase.verifyError(@() multiple.modifier_fun([1;2],17,problem,0),'MATLAB:Feature:CompositeModifier');
        end
    end
end
