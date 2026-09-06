classdef TestPdfToolRuntime < matlab.unittest.TestCase
%TESTPDFTOOLRUNTIME Narrow Linux loader recovery, with real child processes.
    properties
        RunTool
        WorkDir
        OriginalPath
        OriginalLoader
        EnvNames
        EnvValues
    end
    methods (TestMethodSetup)
        function isolate(testCase)
            testCase.assumeTrue(isunix && ~ismac, 'Linux-only loader recovery.');
            testCase.OriginalPath=path;
            testCase.OriginalLoader=getenv('LD_LIBRARY_PATH');
            testCase.EnvNames={'OP_PDF_TEST_CALLS','OP_PDF_TEST_ARGS','OP_PDF_TEST_OUTPUT'};
            testCase.EnvValues=cellfun(@getenv,testCase.EnvNames,'UniformOutput',false);
            testCase.WorkDir=tempname; mkdir(testCase.WorkDir);
            for k=1:numel(testCase.EnvNames)
                setenv(testCase.EnvNames{k},fullfile(testCase.WorkDir,testCase.EnvNames{k}));
            end
            source=fullfile(fileparts(mfilename('fullpath')),'..','..','src');
            addpath(source); old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(source,'private')); testCase.RunTool=@runPdfToolCommand;
            setenv('LD_LIBRARY_PATH','/nonexistent/test-MATLAB-runtime');
        end
    end
    methods (TestMethodTeardown)
        function restore(testCase)
            if isempty(testCase.WorkDir), return; end
            setenv('LD_LIBRARY_PATH',testCase.OriginalLoader);
            for k=1:numel(testCase.EnvNames)
                setenv(testCase.EnvNames{k},testCase.EnvValues{k});
            end
            path(testCase.OriginalPath);
            rmdir(testCase.WorkDir,'s');
        end
    end
    methods (Test)
        function retryPreservesQuotedArgumentsAndParentEnvironment(testCase)
            script=testCase.writeTool([testCase.loaderFailure('GLIBCXX_3.4.32'), ...
                {'printf ''%s\n'' "$@" > "$OP_PDF_TEST_ARGS"', ...
                 'printf ''fresh PDF test payload'' > "$OP_PDF_TEST_OUTPUT"'}]);
            args={fullfile(testCase.WorkDir,'first source ''quote'' $cash ;.pdf'), ...
                fullfile(testCase.WorkDir,'second "quote" $(not_a_command).pdf')};
            command=strjoin([{quote(script)},cellfun(@quote,args,'UniformOutput',false)],' ');
            before=getenv('LD_LIBRARY_PATH'); before_path=getenv('PATH'); before_pwd=pwd;
            [status,~,retried]=testCase.RunTool(command,getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyEqual(status,0); testCase.verifyTrue(retried);
            testCase.verifyEqual(fileread(getenv('OP_PDF_TEST_ARGS')),sprintf('%s\n',args{:}));
            testCase.verifyEqual(fileread(getenv('OP_PDF_TEST_OUTPUT')),'fresh PDF test payload');
            testCase.verifyEqual(testCase.callCount(),2);
            testCase.verifyEqual(getenv('LD_LIBRARY_PATH'),before);
            testCase.verifyEqual(getenv('PATH'),before_path);
            testCase.verifyEqual(pwd,before_pwd);
        end
        function cxxabiLoaderFailureAlsoRetries(testCase)
            script=testCase.writeTool([testCase.loaderFailure('CXXABI_1.3.15'),{'exit 0'}]);
            [status,~,retried]=testCase.RunTool(quote(script),getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyEqual(status,0); testCase.verifyTrue(retried);
            testCase.verifyEqual(testCase.callCount(),2);
        end
        function syntaxAndOrdinaryErrorsDoNotRetry(testCase)
            messages={'Syntax Error: invalid PDF trailer', 'qpdf: command not found', ...
                'version GLIBCXX_3.4.32 not found in a PDF text field'};
            for k=1:numel(messages)
                script=testCase.writeTool({['printf ''%s\n'' ',quote(messages{k}),' >&2'],'exit 65'});
                [status,output,retried]=testCase.RunTool(quote(script),getenv('OP_PDF_TEST_OUTPUT'));
                testCase.verifyEqual(status,65); testCase.verifyFalse(retried);
                testCase.verifySubstring(output,messages{k});
                testCase.verifyEqual(testCase.callCount(),k);
            end
        end
        function missingExecutableDoesNotRetry(testCase)
            [status,~,retried]=testCase.RunTool(quote(fullfile(testCase.WorkDir,'missing tool')), ...
                getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyNotEqual(status,0); testCase.verifyFalse(retried);
            testCase.verifyFalse(isfile(getenv('OP_PDF_TEST_CALLS')));
        end
        function successfulCommandWithDiagnosticTextDoesNotRetry(testCase)
            script=testCase.writeTool({testCase.loaderDiagnostic('GLIBCXX_3.4.32'),'exit 0'});
            [status,~,retried]=testCase.RunTool(quote(script),getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyEqual(status,0); testCase.verifyFalse(retried);
            testCase.verifyEqual(testCase.callCount(),1);
        end
        function absentLoaderSearchPathDoesNotRetry(testCase)
            setenv('LD_LIBRARY_PATH','');
            script=testCase.writeTool({testCase.loaderDiagnostic('GLIBCXX_3.4.32'),'exit 127'});
            [status,~,retried]=testCase.RunTool(quote(script),getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyEqual(status,127); testCase.verifyFalse(retried);
            testCase.verifyEqual(testCase.callCount(),1);
            testCase.verifyEmpty(getenv('LD_LIBRARY_PATH'));
        end
        function failedRetryKeepsBothDiagnosticsAndRemovesOnlyStaging(testCase)
            script=testCase.writeTool([testCase.loaderFailure('GLIBCXX_3.4.32'), ...
                {'printf ''isolated invalid PDF\n'' >&2','exit 65'}]);
            target=fullfile(testCase.WorkDir,'previous summary.pdf'); writeText(target,'old aggregate');
            before=getenv('LD_LIBRARY_PATH');
            [status,output,retried]=testCase.RunTool(quote(script),getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyEqual(status,65); testCase.verifyTrue(retried);
            testCase.verifySubstring(output,'GLIBCXX_3.4.32');
            testCase.verifySubstring(output,'isolated invalid PDF');
            testCase.verifySubstring(output,'Retry without LD_LIBRARY_PATH failed');
            testCase.verifyFalse(isfile(getenv('OP_PDF_TEST_OUTPUT')));
            testCase.verifyEqual(fileread(target),'old aggregate');
            testCase.verifyEqual(testCase.callCount(),2);
            testCase.verifyEqual(getenv('LD_LIBRARY_PATH'),before);
        end
        function successWithoutNewOutputCannotReuseFirstPartial(testCase)
            script=testCase.writeTool([testCase.loaderFailure('GLIBCXX_3.4.32'),{'exit 0'}]);
            [status,~,retried]=testCase.RunTool(quote(script),getenv('OP_PDF_TEST_OUTPUT'));
            testCase.verifyEqual(status,0); testCase.verifyTrue(retried);
            % The merger also requires isfile(temp_output), so an exit-zero
            % child that writes nothing cannot publish the first failed copy.
            testCase.verifyFalse(isfile(getenv('OP_PDF_TEST_OUTPUT')));
        end
        function cannotRemovePartialStopsBeforeRetry(testCase)
            folder=fullfile(testCase.WorkDir,'readonly staging'); mkdir(folder);
            output=fullfile(folder,'partial.pdf'); writeText(output,'partial');
            script=testCase.writeTool({testCase.loaderDiagnostic('GLIBCXX_3.4.32'),'exit 127'});
            [status,diagnostic]=system(['chmod 500 ',quote(folder)]);
            testCase.assertEqual(status,0,diagnostic);
            cleanup=onCleanup(@() system(['chmod 700 ',quote(folder)]));
            [status,output_text,retried]=testCase.RunTool(quote(script),output);
            testCase.verifyNotEqual(status,0); testCase.verifyFalse(retried);
            testCase.verifySubstring(output_text,'Cannot remove the failed temporary PDF');
            testCase.verifyEqual(testCase.callCount(),1);
            testCase.verifyEqual(fileread(output),'partial');
        end
    end
    methods (Access=private)
        function file=writeTool(testCase,body)
            file=fullfile(testCase.WorkDir,'tool space ''quoted'' ; $');
            lines=[{'#!/bin/sh','printf ''call\n'' >> "$OP_PDF_TEST_CALLS"'},body];
            writeText(file,sprintf('%s\n',lines{:}));
            [status,diagnostic]=system(['chmod 700 ',quote(file)]);
            testCase.assertEqual(status,0,diagnostic);
        end
        function count=callCount(~)
            count=numel(strfind(fileread(getenv('OP_PDF_TEST_CALLS')),sprintf('call\n')));
        end
    end
    methods (Static, Access=private)
        function lines=loaderFailure(symbol)
            lines={'if [ -n "${LD_LIBRARY_PATH:-}" ]; then', ...
                'printf ''partial from failed launch'' > "$OP_PDF_TEST_OUTPUT"', ...
                TestPdfToolRuntime.loaderDiagnostic(symbol),'exit 127','fi'};
        end
        function line=loaderDiagnostic(symbol)
            diagnostic=['pdfunite: /opt/MATLAB/sys/os/glnxa64/libstdc++.so.6: version `', ...
                symbol,''' not found (required by /usr/lib/libpoppler.so.134)'];
            line=['printf ''%s\n'' ',quote(diagnostic),' >&2'];
        end
    end
end
function value=quote(value)
    value=char("'"+replace(string(value),"'","'\''")+"'");
end
function writeText(file,text)
    fid=fopen(file,'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s',text);
end
