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
