classdef TestFeatureCustomV2 < matlab.unittest.TestCase
% Custom callbacks use observed immediate predecessors and retained seed reads.
    methods (Test)
        function valueSeedReadAndCallbackProbesAreAllServed(testCase)
            noise=struct('name','noisy','options',struct('noise_mode','deterministic', ...
                'noise_map',@(x)1,'noise_type','absolute','noise_level',2));
            custom=struct('name','custom','options',struct( ...
                'mod_fun',@(x,s,p)p.fun(x)+p.fun(x), ...
                'mod_cub',@(x,s,p)p.cub(x)+p.cub(x), ...
                'mod_ceq',@(x,s,p)p.ceq(x)+p.ceq(x)));
            p=Problem(struct('fun',@(x)x(1)^2,'x0',.5, ...
                'cub',@(x)[x(1)-.25;2*x(1)-.75],'ceq',@(x)x(1)-.5));
            fp=FeaturedProblem(p,Feature({noise,custom}),4,17);
            testCase.verifyEqual(fp.fun(.5),4.5);
            testCase.verifyEqual(fp.cub(.5),[4.5;4.5]);
            testCase.verifyEqual(fp.ceq(.5),4);
            receipt=fp.runtimeReceipt();
            testCase.verifyEqual(receipt.stages{1}.served,struct('fun',3,'cub',3,'ceq',3));
            testCase.verifyEqual(receipt.stages{2}.served,struct('fun',1,'cub',1,'ceq',1));
            testCase.verifyEqual([fp.n_eval_fun,fp.n_eval_cub,fp.n_eval_ceq],[1,1,1]);
            testCase.verifyEqual(fp.fun_hist,.25);
            testCase.verifyEqual(fp.cub_hist,[.25;.25]);
            testCase.verifyEqual(fp.ceq_hist,0);
            fp.evaluateTruth(.5); testCase.verifyEqual(fp.runtimeReceipt(),receipt);
        end
        function customViolationProbeUsesObservedPredecessor(testCase)
            noise=struct('name','noisy','options',struct('noise_mode','deterministic', ...
                'noise_map',@(x)1,'noise_type','absolute','noise_level',2));
            custom=struct('name','custom','options',struct('mod_fun',@(x,s,p)p.maxcv(x)+p.fun(x)));
            p=Problem(struct('fun',@(x)x(1)^2,'x0',.5, ...
                'cub',@(x)[x(1)-.25;2*x(1)-.75],'ceq',@(x)x(1)-.5));
            fp=FeaturedProblem(p,Feature({noise,custom}),4,17);
            testCase.verifyEqual(fp.fun(.5),4.5);
            testCase.verifyEqual(fp.fun_hist,.25);
            testCase.verifyEqual(fp.maxcv_hist,.25);
            receipt=fp.runtimeReceipt();
            testCase.verifyEqual(receipt.stages{1}.served,struct('fun',2,'cub',1,'ceq',1));
            testCase.verifyEqual([fp.n_eval_cub,fp.n_eval_ceq],[0,0]);
        end
        function callbackOutputPoliciesCoverComplexEmptyAndRows(testCase)
            original=warning; cleanup=onCleanup(@()warning(original));
            warning('off','MATLAB:Feature:InvalidCustomObjective');
            p=Problem(struct('fun',@(x)sum(x.^2),'x0',[.5;.5], ...
                'cub',@(x)[x(1)-1;x(2)-1]));
            for value={[],[1,2],[1;2],1+1i,'a',{1}}
                custom=struct('name','custom','options',struct('mod_fun',@(x,s,p)value{1}));
                fp=FeaturedProblem(p,Feature({'noisy',custom}),3,17);
                testCase.verifyTrue(isnan(fp.fun([.5;.5])));
                testCase.verifyEqual(fp.fun_hist,.5);
            end
            for value={NaN,Inf,-Inf,true,single(.5)}
                custom=struct('name','custom','options',struct('mod_fun',@(x,s,p)value{1}));
                fp=FeaturedProblem(p,Feature({'noisy',custom}),3,17);
                testCase.verifyEqual(fp.fun([.5;.5]),double(value{1}));
            end
            custom=struct('name','custom','options',struct('mod_cub',@(x,s,p)[1,2]));
            fp=FeaturedProblem(p,Feature({'noisy',custom}),3,17);
            testCase.verifyEqual(fp.cub([.5;.5]),[1;2]);
            for value={[],[1;2;3],[1,2;3,4],[1;1i],'ab'}
                custom.options.mod_cub=@(x,s,p)value{1};
                fp=FeaturedProblem(p,Feature({'noisy',custom}),3,17);
                testCase.verifyError(@()fp.cub([.5;.5]),'MATLAB:Feature:InvalidCustomOutput');
            end
            % Retained custom kernel skips callbacks on absent channels.
            empty=Problem(struct('fun',@(x)x(1)^2,'x0',.5));
            custom.options.mod_cub=@(x,s,p)error('AB:UnexpectedCallback','Empty callback must be skipped.');
            fp=FeaturedProblem(empty,Feature({'noisy',custom}),3,17);
            testCase.verifyEmpty(fp.cub(.5));
            receipt=fp.runtimeReceipt();
            testCase.verifyEqual(receipt.stages{1}.served.cub,1);
        end
        function constructionShapesAndCallbackErrorsStayExplicit(testCase)
            p=Problem(struct('fun',@(x)sum(x.^2),'x0',[.5;.5]));
            bad={struct('mod_x0',@(s,p)[]), ...
                struct('mod_x0',@(s,p)[1;1i]), ...
                struct('mod_bounds',@(s,p)deal([0,0,0],[1,1,1])), ...
                struct('mod_linear_ub',@(s,p)deal([1,1,1],1)), ...
                struct('mod_linear_eq',@(s,p)deal([1,1],[1;2])), ...
                struct('mod_affine',@(s,p)deal(eye(3),zeros(3,1),eye(3))), ...
                struct('mod_affine',@(s,p)deal(eye(2),0,eye(2)))};
            for k=1:numel(bad)
                spec=Feature({'noisy',struct('name','custom','options',bad{k})});
                testCase.verifyError(@()FeaturedProblem(p,spec,3,17),'MATLAB:Feature:InvalidCustomOutput');
            end
            bad_inverse=struct('name','custom','options',struct('mod_affine',@(s,p)deal(eye(2),zeros(2,1),2*eye(2))));
            testCase.verifyError(@()FeaturedProblem(p,Feature({'noisy',bad_inverse}),3,17), ...
                'MATLAB:Feature:AffineTransformationNotInvertible');
            custom=struct('name','custom','options',struct('mod_bounds',@(s,p)deal([0,0],[1,1]), ...
                'mod_linear_ub',@(s,p)deal([1,1],2)));
            fp=FeaturedProblem(p,Feature({'noisy',custom}),3,17);
            testCase.verifyEqual(fp.xl,[0;0]); testCase.verifyEqual(fp.xu,[1;1]);
            testCase.verifyEqual(fp.aub,[1,1]); testCase.verifyEqual(fp.bub,2);
            custom.options=struct('mod_fun',@(x,s,p)error('AB:UserCallback','User failure stays visible.'));
            fp=FeaturedProblem(p,Feature({'noisy',custom}),3,17);
            testCase.verifyError(@()fp.fun([.5;.5]),'AB:UserCallback');
        end
    end
end
