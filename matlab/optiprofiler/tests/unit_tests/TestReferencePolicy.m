classdef TestReferencePolicy < matlab.unittest.TestCase
% Two public reference kinds, source coordinates, and empirical fallback.
    methods (Test)
        function twoKindsNeedNoMapping(testCase)
            feasible = struct('fun',0,'maxcv',0,'kind','feasible','source','analytic');
            merit = struct('fun',-1,'maxcv',0.25,'kind','merit','source','author');
            p = Problem(struct('fun',@(x)sum(x.^2),'x0',2,'reference',feasible));
            testCase.verifyEqual(p.reference,feasible);
            p = Problem(struct('fun',@(x)sum(x.^2),'x0',2,'reference',merit));
            testCase.verifyEqual(p.reference,merit);
        end
        function builtinsKeepSourceCoordinatesAndNonzeroViolation(testCase)
            record = struct('fun',5,'maxcv',0.25,'kind','merit','source','author','x_ref',[1;2]);
            p = Problem(struct('fun',@(x)sum(x.^2),'x0',[3;4],'reference',record));
            for name = {'permuted','linearly_transformed','noisy+linearly_transformed','linearly_transformed+noisy'}
                q = FeaturedProblem(p,Feature(name{1}),10,17);
                testCase.verifyEqual(q.reference,record);
                testCase.verifyEqual(q.n_eval_fun,0);
            end
        end
        function oldKindsNormalizeWithoutAMappingField(testCase)
            old = {'lower_bound','optimum','best_known','target'};
            canonical = {'feasible','feasible','merit','merit'};
            for k=1:numel(old)
                record = struct('fun',0,'maxcv',0,'kind',old{k},'source','old provider', ...
                    'mapping','objective_violation/1');
                p = Problem(struct('fun',@(x)sum(x.^2),'x0',2,'reference',record));
                testCase.verifyEqual(p.reference.kind,canonical{k});
                testCase.verifyFalse(isfield(p.reference,'mapping'));
            end
        end
        function unsupportedStagesDoNotRestoreAReference(testCase)
            record = struct('fun',0,'maxcv',0,'kind','feasible','source','author','x_ref',[0;0]);
            p = Problem(struct('fun',@(x)sum(x.^2),'x0',[3;4],'reference',record));
            for name = {'custom','noisy+custom','custom+permuted','quantized','quantized+noisy'}
                q = FeaturedProblem(p,Feature(name{1}),10,17);
                testCase.verifyEmpty(q.reference);
            end
        end
        function publicReportExplainsFirstFallback(testCase)
            record = struct('fun',-1,'maxcv',0,'kind','feasible','source','author','x_ref',0);
            p = Problem(struct('fun',@(x)sum(x.^2),'x0',2,'reference',record));
            names = {'custom','quantized+custom','custom+quantized'};
            reasons = {'custom_feature','quantized_ground_truth','custom_feature'};
            for k=1:numel(names)
                file = [tempname,'.json'];
                guard = onCleanup(@() delete(file)); %#ok<NASGU>
                benchmark({@TestReferencePolicy.stay,@TestReferencePolicy.zero}, ...
                    struct('problem',p,'feature_name',names{k},'score_only',true,'n_jobs',1, ...
                    'n_runs',1,'run_plain',false,'silent',true,'report_path',file));
                report = jsondecode(fileread(file));
                testCase.verifyEqual(report.problems.reference_fallback_reasons,{reasons{k}});
            end
        end
    end
    methods (Static)
        function x = stay(fun,x0)
            fun(x0); x=x0;
        end
        function x = zero(fun,x0)
            x=zeros(size(x0)); fun(x);
        end
    end
end
