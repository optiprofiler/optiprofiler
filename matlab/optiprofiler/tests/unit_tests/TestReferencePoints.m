classdef TestReferencePoints < matlab.unittest.TestCase
% Optional reference witnesses survive persistence and coordinate changes.
    methods (Static)
        function record = reference()
            record = struct('fun',5,'maxcv',0,'kind','merit', ...
                'source','analytic','x_ref',[1;2]);
        end
        function [A,b,inverse] = affine(~,~)
            A = diag([2,4]); b = [-1;1]; inverse = diag([0.5,0.25]);
        end
        function write(file, text)
            fid = fopen(file,'w'); assert(fid>=0);
            guard = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fprintf(fid,'%s',text);
        end
        function x = stay(fun,x0)
            fun(x0); x = x0;
        end
        function x = zero(fun,x0)
            x = zeros(size(x0)); fun(x);
        end
        function value = bytes(file)
            fid=fopen(file,'rb'); assert(fid>=0);
            guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
            value=fread(fid,Inf,'*uint8');
        end
        function cleanup(root,registry,old_path,old_cwd)
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',registry);
            cd(old_cwd);
            path(old_path);
            if isfolder(root), rmdir(root,'s'); end
        end
    end
    methods (Test)
        function pointIsOptionalAndPersisted(testCase)
            record = TestReferencePoints.reference();
            problem = Problem(struct('fun',@(x) sum(x.^2),'x0',[3;4],'reference',record));
            testCase.verifyEqual(problem.reference.x_ref,[1;2]);
            file = [tempname,'.mat'];
            cleanup = onCleanup(@() delete(file)); %#ok<NASGU>
            save(file,'problem');
            restored = load(file,'problem');
            testCase.verifyEqual(restored.problem.reference,record);
        end
        function pointStaysInOriginalCoordinates(testCase)
            record = TestReferencePoints.reference();
            problem = Problem(struct('fun',@(x) sum(x.^2),'x0',[3;4],'reference',record));
            features = {Feature('linearly_transformed'), Feature('linearly_transformed+noisy')};
            for k = 1:numel(features)
                featured = FeaturedProblem(problem,features{k},10,17);
                testCase.verifyEqual(featured.reference.x_ref,[1;2]);
                testCase.verifyEqual(featured.reference.fun,5);
                testCase.verifyEqual(featured.n_eval_fun,0);
            end
            testCase.verifyEqual(problem.reference.x_ref,[1;2]);
        end
        function pointDimensionMustMatch(testCase)
            record = TestReferencePoints.reference(); record.x_ref = 1;
            testCase.verifyError(@() Problem(struct('fun',@(x) sum(x.^2),'x0',[3;4], ...
                'reference',record)), 'MATLAB:Problem:reference_x_ref_Dimension');
        end
        function publicArchiveReloadRetainsWitness(testCase)
            root = tempname; mkdir(root);
            old_cwd = pwd;
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            old_path = path;
            guard = onCleanup(@() TestReferencePoints.cleanup(root,registry,old_path,old_cwd)); %#ok<NASGU>
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',fullfile(root,'registry.mat'));
            library = fullfile(root,'library'); mkdir(library);
            TestReferencePoints.write(fullfile(library,'refpoints_select.m'), sprintf( ...
                'function names=refpoints_select(~)\nnames={''P''};\nend\n'));
            TestReferencePoints.write(fullfile(library,'refpoints_load.m'), sprintf([ ...
                'function problem=refpoints_load(name)\n' ...
                'r=struct(''fun'',0,''maxcv'',0,''kind'',''optimum'',''source'',''analytic'',''mapping'',''objective_violation/1'',''x_ref'',0);\n' ...
                'problem=Problem(struct(''fun'',@(x)sum(x.^2),''x0'',2,''name'',name,''reference'',r));\nend\n']));
            registerProblemLibrary(struct('name','refpoints','root',library, ...
                'select_function','refpoints_select','load_function','refpoints_load'));
            options = struct('plibs',{{'refpoints'}},'n_jobs',1,'n_runs',1,'max_eval_factor',2, ...
                'max_tol_order',1,'run_plain',false,'draw_hist_plots','none','silent',true, ...
                'savepath',root,'benchmark_id','witness','report_path',fullfile(root,'report.json'));
            for name={'plain','custom+noisy'}
                options.feature_name=name{1};
                options.benchmark_id=['witness-',strrep(name{1},'+','-')];
                options.report_path=fullfile(root,[options.benchmark_id,'.json']);
                original = benchmark({@TestReferencePoints.stay,@TestReferencePoints.zero},options);
                files = dir(fullfile(root,options.benchmark_id,'**','data_for_loading.mat'));
                testCase.assertNumElements(files,1);
                archive=fullfile(files(1).folder,files(1).name);
                original_bytes=TestReferencePoints.bytes(archive);
                stored = load(archive,'results_plibs');
                if strcmp(name{1},'plain')
                    testCase.verifyEqual(stored.results_plibs{1}.reference_records{1}.x_ref,0);
                    raw = fileread(options.report_path);
                    testCase.verifyNotEmpty(regexp(raw,'"x_ref"\s*:\s*\[\s*0\s*\]','once'));
                else
                    testCase.verifyEmpty(stored.results_plibs{1}.reference_records{1});
                    testCase.verifyEqual(stored.results_plibs{1}.reference_fallback_reasons{1},'custom_feature');
                end
                reload=fullfile(root,['reload-',options.benchmark_id,'.json']);
                % Reload source is cwd/benchmark_id, independently of output.
                cd(root);
                loaded = benchmark(struct('load','latest','savepath',root,'benchmark_id',options.benchmark_id, ...
                    'score_only',true,'silent',true,'max_tol_order',1,'report_path',reload));
                cd(old_cwd);
                testCase.verifyEqual(loaded,original);
                testCase.verifyEqual(TestReferencePoints.bytes(archive),original_bytes);
                if ~strcmp(name{1},'plain')
                    report=jsondecode(fileread(reload));
                    testCase.verifyEqual(report.problems.reference_fallback_reasons,{'custom_feature'});
                end
            end
        end
    end
end
