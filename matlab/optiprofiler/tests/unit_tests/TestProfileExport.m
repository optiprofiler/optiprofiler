classdef TestProfileExport < matlab.unittest.TestCase
% Public benchmark must export tiled summaries and append tolerance pages.
    methods (Test)
        function publicBenchmarkExportsTiledSummaries(testCase)
            folder=fileparts(mfilename('fullpath'));
            root=tempname; mkdir(root);
            old=pwd; oldpath=path; oldrng=rng;
            registry=getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            testCase.addTeardown(@() rmdir(root,'s'));
            testCase.addTeardown(@() rng(oldrng));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',registry));
            testCase.addTeardown(@() path(oldpath));
            testCase.addTeardown(@() cd(old));
            addpath(fullfile(folder,'..','..','src'));
            fixture=fullfile(folder,'..','fixtures','profilebands'); addpath(fixture);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',fullfile(root,'registry.mat'));
            registerProblemLibrary(struct('name','export_fixture','root',fixture, ...
                'select_function','band_fixture_select','load_function','band_fixture_load'));
            cd(root);
            options=struct('plibs',{{'export_fixture'}},'ptype','u','mindim',1,'maxdim',10, ...
                'n_jobs',1,'n_runs',1,'max_eval_factor',2,'max_tol_order',2, ...
                'draw_hist_plots','none','silent',true,'benchmark_id','exported', ...
                'solver_names',{{'stay','trial'}});
            benchmark({@stay,@trial},options);
            pages=dir(fullfile(root,'exported','**','data_out_*.pdf'));
            testCase.assertEqual(numel(pages),2,'Both tolerance pages must be exported.');
            combined=dir(fullfile(root,'exported','**','data_out.pdf'));
            testCase.assertEqual(numel(combined),1);
            testCase.verifyGreaterThan(combined.bytes,max([pages.bytes]));
            summaries=dir(fullfile(root,'exported','**','summary_*.pdf'));
            testCase.assertNotEmpty(summaries,'The nested combined layout must be exported.');
            for file=[pages(:);combined(:);summaries(:)]'
                fid=fopen(fullfile(file.folder,file.name),'r');
                header=char(fread(fid,5,'uint8')'); fclose(fid);
                testCase.verifyEqual(header,'%PDF-');
            end
        end
    end
end

function x=stay(fun,x0)
    fun(x0); x=x0;
end
function x=trial(fun,x0)
    fun(x0); x=2*x0; fun(x);
end
