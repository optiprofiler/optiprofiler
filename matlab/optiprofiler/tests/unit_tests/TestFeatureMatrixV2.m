classdef TestFeatureMatrixV2 < matlab.unittest.TestCase
% Bounded public matrix with literal oracles, not regenerated golden values.
    methods (Test)
        function everyOrderedBuiltInPairAcrossProblemTypes(testCase)
            entries = TestFeatureMatrixV2.entries();
            types = {'u','b','l','n'};
            original_rng = rng;
            cases = 0;
            for kind = 1:numel(types)
                root = struct('fun',@(x) 2,'x0',[0;0]);
                if kind >= 2
                    root.xl = [-1;-2]; root.xu = [1;2];
                end
                if kind >= 3
                    root.aub = [1,1;-1,2]; root.bub = [4;4];
                    root.aeq = [0,0]; root.beq = 0;
                end
                if kind == 4
                    root.cub = @(x) [-1;-2]; root.ceq = @(x) 0;
                end
                problem = Problem(root);
                testCase.verifyEqual(problem.ptype,types{kind});
                for first = 1:10
                    for second = 1:10
                        spec = Feature({entries{first},entries{second}});
                        fp = FeaturedProblem(problem,spec,3,17);
                        x = fp.x0;
                        [fun,cub,ceq] = TestFeatureMatrixV2.literalPair(first,second,kind==4);
                        testCase.verifyEqual(fp.fun(x),fun);
                        testCase.verifyEqual(fp.cub(x),cub);
                        testCase.verifyEqual(fp.ceq(x),ceq);
                        testCase.verifyEqual(fp.fun_init,2);
                        testCase.verifyEqual(fp.fun_hist,2);
                        testCase.verifyEqual(fp.maxcv_init,0);
                        testCase.verifyEqual(fp.maxcv_hist,0);
                        if kind == 4
                            testCase.verifyEqual(fp.cub_hist,[-1;-2]);
                            testCase.verifyEqual(fp.ceq_hist,0);
                            testCase.verifyEqual([fp.m_nonlinear_ub,fp.m_nonlinear_eq],[2,1]);
                        else
                            testCase.verifySize(fp.cub_hist,[0,1]);
                            testCase.verifySize(fp.ceq_hist,[0,1]);
                            testCase.verifyEqual([fp.m_nonlinear_ub,fp.m_nonlinear_eq],[0,0]);
                        end
                        testCase.verifyEqual([fp.n_eval_fun,fp.n_eval_cub,fp.n_eval_ceq],[1,1,1]);
                        testCase.verifySize(fp.x0,[2,1]);
                        testCase.verifySize(fp.xl,[2,1]); testCase.verifySize(fp.xu,[2,1]);
                        testCase.verifyEqual(size(fp.aub,2),2); testCase.verifyEqual(size(fp.aeq,2),2);
                        receipt = fp.runtimeReceipt();
                        fp.evaluateTruth(x); fp.maxcv(x,true); fp.toOriginalCoordinates(x);
                        testCase.verifyEqual(fp.runtimeReceipt(),receipt);
                        testCase.verifyEqual(numel(receipt.stages),2);
                        testCase.verifyEqual(fp.execution_strategy,'composed-views');
                        cases = cases+1;
                    end
                end
            end
            testCase.verifyEqual(cases,400);
            testCase.verifyEqual(rng,original_rng);
        end
        function mixedThirtyTwoAndSixtyFourStageChainsUseSameContracts(testCase)
            entries = TestFeatureMatrixV2.entries();
            entries{6}.options.nan_rate = 0;
            p = Problem(struct('fun',@(x) 2,'x0',[0;0]));
            for count = [32,64]
                indices = mod(0:count-1,10)+1;
                fp = FeaturedProblem(p,Feature(entries(indices)),3,17);
                % Only literal +1/4 noise and +1/2 custom stages alter this
                % constant unconstrained oracle; meshes/maps preserve it.
                expected = 2+sum(indices==2)/4+sum(indices==10)/2;
                testCase.verifyEqual(fp.fun(fp.x0),expected);
                testCase.verifyEqual(fp.fun_hist,2);
                testCase.verifyEqual(fp.fun_init,2);
                testCase.verifyEqual(fp.maxcv_hist,0);
                receipt = fp.runtimeReceipt();
                testCase.verifyEqual(numel(receipt.stages),count);
                fp.evaluateTruth(fp.x0);
                testCase.verifyEqual(fp.runtimeReceipt(),receipt);
            end
        end
        function activeBuiltInAffineAndMeshOrderHasHandComputedTruth(testCase)
            p = Problem(struct('fun',@(x)[1,3]*x,'x0',[1;1], ...
                'xl',[0;0],'xu',[2;3],'aub',[1,2],'bub',5, ...
                'aeq',[1,-1],'beq',0,'cub',@(x)sum(x)-1,'ceq',@(x)x(1)-x(2)));
            affine = struct('name','linearly_transformed','options', ...
                struct('rotated',false,'condition_factor',4));
            x = [1.25;.375];
            for truth = [false,true]
                quantized = struct('name','quantized','options',struct('mesh_size',.5,'ground_truth',truth));
                a = FeaturedProblem(p,Feature({affine,quantized}),3,17);
                b = FeaturedProblem(p,Feature({quantized,affine}),3,17);
                testCase.verifyEqual(a.x0,[2;.5]);
                testCase.verifyEqual(a.xl,[0;0]); testCase.verifyEqual(a.xu,[4;1.5]);
                testCase.verifyEqual(a.aub,[.5,4]); testCase.verifyEqual(a.bub,5);
                testCase.verifyEqual(a.aeq,[.5,-2]); testCase.verifyEqual(a.beq,0);
                testCase.verifyEqual(a.fun(x),3.75); testCase.verifyEqual(b.fun(x),3.5);
                testCase.verifyEqual(a.toOriginalCoordinates(x),[.625;.75]);
                testCase.verifyEqual(b.toOriginalCoordinates(x),[.625;.75]);
                if truth
                    testCase.verifyEqual(a.fun_hist,3.75); testCase.verifyEqual(b.fun_hist,3.5);
                    testCase.verifyEqual(a.maxcv_hist,.75); testCase.verifyEqual(b.maxcv_hist,.5);
                else
                    testCase.verifyEqual(a.fun_hist,2.875); testCase.verifyEqual(b.fun_hist,2.875);
                    testCase.verifyEqual(a.maxcv_hist,.375); testCase.verifyEqual(b.maxcv_hist,.375);
                end
            end
        end
    end
    methods (Static)
        function entries = entries()
            entries = {struct('name','perturbed_x0','options',struct('perturbation_level',1e-3)), ...
                struct('name','noisy','options',struct('noise_mode','deterministic','noise_map',@(x)1, ...
                    'noise_level',.25,'noise_type','absolute')), ...
                struct('name','truncated','options',struct('significant_digits',6)), ...
                struct('name','permuted'), ...
                struct('name','linearly_transformed','options',struct('rotated',false,'condition_factor',4)), ...
                struct('name','random_nan','options',struct('nan_rate',1)), ...
                struct('name','unrelaxable_constraints','options',struct('unrelaxable_bounds',true, ...
                    'unrelaxable_linear_constraints',true,'unrelaxable_nonlinear_constraints',true)), ...
                struct('name','nonquantifiable_constraints'), ...
                struct('name','quantized','options',struct('mesh_size',.125,'ground_truth',true)), ...
                struct('name','custom','options',struct('mod_fun',@(x,s,p)p.fun(x)+.5, ...
                    'mod_cub',@(x,s,p)2*p.cub(x),'mod_ceq',@(x,s,p)p.ceq(x)+.125))};
        end
        function [fun,cub,ceq] = literalPair(first,second,nonlinear)
            fun = 2;
            cub = zeros(0,1); ceq = zeros(0,1);
            if nonlinear
                cub = [-1;-2]; ceq = 0;
            end
            % Hand-defined transformations of these fixed literal values.
            % No Feature, runtime, RNG, reference or candidate kernel calls.
            for index = [first,second]
                switch index
                    case 2
                        fun = fun+.25; cub = cub+.25; ceq = ceq+.25;
                    case 6
                        fun = NaN; cub(:) = NaN; ceq(:) = NaN;
                    case 7
                        if any(cub>0) || any(abs(ceq)>0)
                            fun = Inf;
                        end
                    case 8
                        cub(cub<=0) = 0; cub(cub>0) = 1;
                        ceq(abs(ceq)<=1e-6) = 0; ceq(abs(ceq)>1e-6) = 1;
                    case 10
                        fun = fun+.5; cub = 2*cub; ceq = ceq+.125;
                end
            end
        end
    end
end
