classdef TestFeatureCompositionV2 < matlab.unittest.TestCase
% Public numeric oracle/reference behavior of ordered problem views.
    methods (Test)
        function repeatedNoiseStagesHaveObservedAndReferenceOutputs(testCase)
            feature = Feature({struct('name','noisy','options',struct( ...
                'noise_mode','deterministic','noise_map',@(x) 1, ...
                'noise_type','absolute','noise_level',2)), ...
                struct('name','noisy','options',struct('noise_mode','deterministic', ...
                'noise_map',@(x) 1,'noise_type','absolute','noise_level',3))});
            problem = Problem(struct('fun',@(x) x(1)^2,'x0',.5));
            featured = FeaturedProblem(problem,feature,4,17);
            testCase.verifyEqual(featured.fun_init,.25);
            testCase.verifyEqual(featured.fun(.5),5.25);
            testCase.verifyEqual(featured.fun_hist,.25);
            [truth,violation] = featured.evaluateTruth(.5);
            testCase.verifyEqual(truth,.25);
            testCase.verifyEqual(violation,0);
            testCase.verifyEqual(featured.n_eval_fun,1);
            testCase.verifyEqual(featured.x0,.5);
        end
        function quantizationHasStageLocalReferenceAndLazyReads(testCase)
            noise=struct('name','noisy','options',struct('noise_mode','deterministic', ...
                'noise_map',@(x) 1,'noise_type','absolute','noise_level',2));
            quantized=struct('name','quantized','options',struct('mesh_size',.5,'ground_truth',true));
            p=Problem(struct('fun',@(x)x(1)^2,'x0',.6,'xl',0,'xu',.55, ...
                'aub',1,'bub',.58,'cub',@(x)[x(1)-.2;2*x(1)-1],'ceq',@(x)x(1)-.5));
            fp=FeaturedProblem(p,Feature({noise,quantized}),4,17);
            testCase.verifyEqual(fp.fun_init,.25);
            testCase.verifyEqual(fp.fun(.6),2.25);
            testCase.verifyEqual(fp.fun_hist,.25);
            testCase.verifyEqual(fp.cub(.6),[2.3;2]);
            testCase.verifyEqual(fp.cub_hist,[.3;0]);
            testCase.verifyEqual(fp.n_eval_cub,1);
            [cv,bounds,linear,nonlinear]=fp.maxcv(.6,true);
            testCase.verifyEqual(cv,.3);
            testCase.verifyEqual(bounds,.6-.55);
            testCase.verifyEqual(linear,.6-.58);
            testCase.verifyEqual(nonlinear,.3);
            before=fp.runtimeReceipt();
            fp.evaluateTruth(.6); fp.maxcv(.6);
            testCase.verifyEqual(fp.runtimeReceipt(),before);
            testCase.verifyEqual(before.stages{1}.served.fun,1);
            testCase.verifyEqual(before.stages{1}.served.cub,1);
            testCase.verifyEqual(fp.toOriginalCoordinates(.6),.6);
            quantized.options.ground_truth=false;
            other=FeaturedProblem(p,Feature({noise,quantized}),4,17);
            testCase.verifyEqual(other.fun_init,.6^2);
            testCase.verifyEqual(other.fun(.6),2.25);
            testCase.verifyEqual(other.fun_hist,.6^2);
        end
        function affineOrderTransportsBoundsLinearsAndOutput(testCase)
            first=struct('name','custom','options',struct('mod_affine',@TestFeatureCompositionV2.firstAffine));
            second=struct('name','custom','options',struct('mod_affine',@TestFeatureCompositionV2.secondAffine));
            p=Problem(struct('fun',@(x)[1,10]*x,'x0',[1;1], ...
                'xl',[0;0],'xu',[2;3],'aub',[1,1],'bub',3,'ceq',@(x)sum(x)-2));
            fp=FeaturedProblem(p,Feature({first,second}),4,17);
            testCase.verifyEqual(fp.x0,[0;1]);
            testCase.verifyEqual(fp.xl,[-Inf;-Inf]);
            testCase.verifyEqual(fp.xu,[Inf;Inf]);
            testCase.verifyEqual(fp.aub,[1,1;0,1;-1,-1;0,-1;2,1]);
            testCase.verifyEqual(fp.bub,[1.5;2;-.5;1;2]);
            testCase.verifyEqual(fp.fun([.5;.25]),18);
            testCase.verifyEqual(fp.toOriginalCoordinates([.5;.25]),[.5;1.75]);
            [truth,cv]=fp.evaluateTruth([.5;.25]);
            testCase.verifyEqual(truth,18);
            testCase.verifyEqual(cv,.25);
            testCase.verifyEqual(fp.m_nonlinear_eq,1);
            testCase.verifyError(@()fp.grad([.5;.25]),'MATLAB:FeaturedProblem:UnsupportedCompositeDerivative');
        end
        function unrelaxableStructuralChecksAreLazyObservedQueries(testCase)
            p=Problem(struct('fun',@(x)x(1)^2,'x0',.25,'xl',0,'xu',.5, ...
                'cub',@(x)[x(1)-.4;x(1)-.3],'ceq',@(x)x(1)-.25));
            noise=struct('name','noisy','options',struct('noise_level',0));
            fp=FeaturedProblem(p,Feature({noise,'unrelaxable_constraints'}),4,17);
            testCase.verifyEqual(fp.fun(.75),Inf);
            receipt=fp.runtimeReceipt();
            testCase.verifyEqual(receipt.stages{1}.served,struct('fun',1,'cub',0,'ceq',0));
            testCase.verifyEqual([fp.n_eval_cub,fp.n_eval_ceq],[0,0]);
            nonlinear=struct('name','unrelaxable_constraints','options', ...
                struct('unrelaxable_bounds',false,'unrelaxable_nonlinear_constraints',true));
            other=FeaturedProblem(p,Feature({noise,nonlinear}),4,17);
            testCase.verifyEqual(other.fun(.375),Inf);
            receipt=other.runtimeReceipt();
            testCase.verifyEqual(receipt.stages{1}.served,struct('fun',1,'cub',1,'ceq',1));
            testCase.verifyEqual([other.n_eval_cub,other.n_eval_ceq],[0,0]);
        end
        function customShapeFailuresStayInObservedChannel(testCase)
            p=Problem(struct('fun',@(x)sum(x.^2),'x0',[1;1], ...
                'cub',@(x)[x(1)-2;x(2)-2],'ceq',@(x)sum(x)-2));
            custom=struct('name','custom','options',struct('mod_fun',@(x,s,p)[1;2]));
            fp=FeaturedProblem(p,Feature({'noisy',custom}),4,17);
            warning_state=warning; cleanup=onCleanup(@()warning(warning_state));
            warning('off','MATLAB:Feature:InvalidCustomObjective');
            testCase.verifyTrue(isnan(fp.fun([1;1])));
            testCase.verifyEqual(fp.fun_hist,2);
            testCase.verifyEqual(fp.fun_init,2);
            custom.options=struct('mod_cub',@(x,s,p)[1,2;3,4]);
            fp=FeaturedProblem(p,Feature({'noisy',custom}),4,17);
            testCase.verifyError(@()fp.cub([1;1]),'MATLAB:Feature:InvalidCustomOutput');
            custom.options=struct('mod_x0',@(s,p)[1;2;3]);
            testCase.verifyError(@()FeaturedProblem(p,Feature({'noisy',custom}),4,17), ...
                'MATLAB:Feature:InvalidCustomOutput');
        end
        function emptyChannelQueriesHaveRecordedQueryBudgets(testCase)
            p=Problem(struct('fun',@(x)x(1)^2,'x0',1,'cub',@(x)[],'ceq',@(x)[]));
            fp=FeaturedProblem(p,Feature('noisy+truncated'),2,17);
            fp.cub(1); fp.ceq(1);
            testCase.verifyEqual([fp.n_eval_cub,fp.n_eval_ceq],[1,1]);
            fp.cub(2,false); fp.ceq(2,false);
            testCase.verifyEqual([fp.n_eval_cub,fp.n_eval_ceq],[1,1]);
            fp.cub(3); fp.ceq(3);
            testCase.verifyEqual([fp.n_eval_cub,fp.n_eval_ceq],[2,2]);
            before=fp.runtimeReceipt(); fp.cub(4); fp.ceq(4);
            testCase.verifyEqual(fp.runtimeReceipt(),before);
            testCase.verifyError(@()fp.cub(5),'MATLAB:FeaturedProblem:cubExceedTerminationEval');
            testCase.verifyError(@()fp.ceq(5),'MATLAB:FeaturedProblem:ceqExceedTerminationEval');
        end
        function allBuiltInKindsWorkInBothOrders(testCase)
            entries={struct('name','perturbed_x0','options',struct('perturbation_level',0)), ...
                struct('name','noisy','options',struct('noise_level',0)), 'truncated','permuted', ...
                struct('name','linearly_transformed','options',struct('rotated',false,'condition_factor',0)), ...
                struct('name','random_nan','options',struct('nan_rate',0)), ...
                'unrelaxable_constraints','nonquantifiable_constraints','quantized','custom'};
            noise=struct('name','noisy','options',struct('noise_level',0));
            p=Problem(struct('fun',@(x)sum(x.^2),'x0',[1;1],'xl',[0;0],'xu',[2;2], ...
                'cub',@(x)sum(x)-3,'ceq',@(x)sum(x)-2));
            original_rng=rng;
            for k=1:numel(entries)
                for reverse=[false,true]
                    stages={entries{k},noise}; if reverse, stages=fliplr(stages); end
                    fp=FeaturedProblem(p,Feature(stages),3,17);
                    testCase.verifyEqual(fp.fun([1;1]),2);
                    fp.cub([1;1]); fp.ceq([1;1]);
                    testCase.verifyEqual(fp.fun_hist,2);
                    testCase.verifyEqual(fp.cub_hist,-1);
                    testCase.verifyEqual(fp.ceq_hist,0);
                    testCase.verifyEqual(fp.maxcv_init,0);
                    receipt=fp.runtimeReceipt();
                    testCase.verifyEqual(numel(receipt.stages),2);
                end
            end
            testCase.verifyEqual(rng,original_rng);
        end
        function orderedNoiseAndRepeatedMeshesKeepLocalTruth(testCase)
            p=Problem(struct('fun',@(x)x(1)^2,'x0',.375));
            noise=struct('name','noisy','options',struct('noise_mode','deterministic', ...
                'noise_map',@(x)x(1),'noise_type','absolute','noise_level',2));
            quantized=struct('name','quantized','options',struct('mesh_size',.5));
            a=FeaturedProblem(p,Feature({noise,quantized}),3,17);
            b=FeaturedProblem(p,Feature({quantized,noise}),3,17);
            testCase.verifyEqual(a.fun(.375),1.25);
            testCase.verifyEqual(b.fun(.375),1);
            testCase.verifyEqual(a.fun_hist,.25); testCase.verifyEqual(b.fun_hist,.25);
            flags=logical([0,0;0,1;1,0;1,1]);
            forward_truth=[9/64,0,1/4,0]; reverse_truth=[9/64,1/4,0,1];
            for k=1:4
                first=struct('name','quantized','options',struct('mesh_size',.25,'ground_truth',flags(k,1)));
                second=struct('name','quantized','options',struct('mesh_size',1,'ground_truth',flags(k,2)));
                a=FeaturedProblem(p,Feature({first,second}),3,17);
                testCase.verifyEqual(a.fun(.375),0);
                testCase.verifyEqual(a.fun_hist,forward_truth(k));
                first.options.mesh_size=1; second.options.mesh_size=.25;
                b=FeaturedProblem(p,Feature({first,second}),3,17);
                testCase.verifyEqual(b.fun(.375),1);
                testCase.verifyEqual(b.fun_hist,reverse_truth(k));
                testCase.verifyEqual(b.toOriginalCoordinates(.375),.375);
            end
        end
        function unrelaxableUsesTransportedConstraintCategories(testCase)
            p=Problem(struct('fun',@(x)sum(x.^2),'x0',[.5;.5], ...
                'xl',[0;0],'xu',[1;1]));
            affine=struct('name','custom','options',struct('mod_affine',@TestFeatureCompositionV2.secondAffine));
            bounds=FeaturedProblem(p,Feature({affine,'unrelaxable_constraints'}),3,17);
            testCase.verifyEqual(bounds.fun([3;0]),4);
            testCase.verifyEqual(bounds.maxcv_hist,1);
            linear=struct('name','unrelaxable_constraints','options', ...
                struct('unrelaxable_bounds',false,'unrelaxable_linear_constraints',true));
            fp=FeaturedProblem(p,Feature({affine,linear}),3,17);
            testCase.verifyEqual(fp.fun([3;0]),Inf);
            testCase.verifyEqual(fp.fun_hist,4);
            testCase.verifyEqual(fp.m_linear_ub,4);
        end
    end
    methods(Static)
        function [A,b,inverse]=firstAffine(~,~)
            A=[2,0;0,-1]; b=[1;2]; inverse=[.5,0;0,-1];
        end
        function [A,b,inverse]=secondAffine(~,~)
            A=[1,1;0,1]; b=[-1;0]; inverse=[1,-1;0,1];
        end
    end
end
