classdef TestReportTextSafety < matlab.unittest.TestCase
%TESTREPORTTEXTSAFETY Report provenance is bounded honestly and keeps no credentials.
%   Mirrors Python's test_report_text_safety.py.
    properties
        OldCwd
        OldPath
        WorkDir
    end
    methods (TestMethodSetup)
        function setupState(testCase)
            testCase.OldCwd = pwd;
            testCase.OldPath = path;
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            cd(testCase.WorkDir);
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.OldCwd);
            path(testCase.OldPath);
            if isfolder(testCase.WorkDir), rmdir(testCase.WorkDir, 's'); end
        end
    end
    methods (Static)
        function x = finish(fun, x0, varargin)
            x = ones(size(x0));
            fun(x);
        end
        function x = stay(fun, x0, varargin)
            fun(x0);
            x = x0;
        end
        function replaceFile(source, target)
            [ok, message] = movefile(source, target, 'f');
            assert(ok, message);
        end
        function value = encoded(value)
            value = jsondecode(optiprofiler_internal.EvalReport.encodeMetadata(value));
        end
        function source = reportedSource(work_dir, source)
            reference = struct('fun', 0, 'maxcv', 0, 'kind', 'feasible', 'source', source);
            options = struct('score_only', true, 'silent', true, 'n_jobs', 1, 'max_eval_factor', 2);
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, 'reference', reference));
            options.report_path = [tempname(work_dir), '.json'];
            benchmark({@TestReportTextSafety.finish, @TestReportTextSafety.stay}, options);
            report = jsondecode(fileread(options.report_path));
            problem = report.problems;
            if iscell(problem), problem = problem{1}; else, problem = problem(1); end
            records = problem.reference_records;
            if iscell(records), record = records{1}; else, record = records(1); end
            source = record.source;
        end
    end
    methods (Test)
        function referenceSourceKeepsLocationWithoutSecrets(testCase)
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                'https://user:pw@doi.org/10.1080/10556788.2025.2490640?utm=mail#top');
            testCase.verifyEqual(reported, ...
                'https://doi.org/10.1080/10556788.2025.2490640 [url parameters removed]');
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                'https://doi.org/10.1080/10556788.2025.2490640');
            testCase.verifyEqual(reported, 'https://doi.org/10.1080/10556788.2025.2490640');
            % Same text as Python's provenance_text for several URLs, a port,
            % letter case and a URL without a host.
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                'Data from HTTP://User@Example.ORG:8080/a/b?x=1 and https://doi.org/10.1/x');
            testCase.verifyEqual(reported, ...
                'Data from http://example.org:8080/a/b [url parameters removed] and https://doi.org/10.1/x');
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, 'https://@/x');
            testCase.verifyEqual(reported, '[redacted_url]');
        end
        function urlOfAnySchemeLosesItsCredentials(testCase)
            canary = 'CANARYpw';
            for scheme = {'https', 'HTTPS', 'Http', 'git+https', 'GIT+SSH', 'git+ssh', 'ssh', 'ftp', 's3', 'custom.scheme-1'}
                url = sprintf('%s://user:%s@example.org/repo.git', scheme{1}, canary);
                testCase.verifyEqual(TestReportTextSafety.encoded(url), '[redacted_url]');
                testCase.verifyEqual(TestReportTextSafety.encoded(string(url)), '[redacted_url]');
                testCase.verifyEqual(TestReportTextSafety.encoded(['built from ', url, ' (pinned)']), ...
                    'built from [redacted_url] (pinned)');
            end
            for scheme = {'HTTPS', 'git+https', 'GIT+SSH', 'ssh', 'custom.scheme-1'}
                url = sprintf('%s://user:%s@example.org/repo.git', scheme{1}, canary);
                expected = [lower(scheme{1}), '://example.org/repo.git [url parameters removed]'];
                reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ['Data from ', url, ' and more']);
                testCase.verifyEqual(reported, ['Data from ', expected, ' and more']);
            end
        end
        function urlGluedToAnotherKeepsNoUserInformation(testCase)
            % The same list as in test_report_text_safety.py;
            % tools/tests/test_report_redaction_contract.py compares them.
            glued = {'https://user:A@example.org/x,https://u2:CANARYpw@example.org/y', ...
                'https://example.org/x,https://u2:CANARYpw@example.org/y', ...
                'https://user:A@example.org,https://u2:CANARYpw@example.org/y', ...
                'https://example.org,https://u2:CANARYpw@example.org/y', ...
                'https://example.org/x;git+ssh://u2:CANARYpw@example.org/y', ...
                'HTTPS://example.org/x,SSH://u2:CANARYpw@example.org/y'};
            for k = 1:numel(glued)
                text = glued{k};
                testCase.verifyEqual(optiprofiler_internal.EvalReport.credentialFreeUrls(text), '[redacted_url]');
                testCase.verifyEqual(optiprofiler_internal.EvalReport.credentialFreeUrls(['Data from ', text, ' and more']), ...
                    'Data from [redacted_url] and more');
                testCase.verifyEqual(TestReportTextSafety.encoded(text), '[redacted_url]');
            end
            testCase.verifyEqual(TestReportTextSafety.reportedSource(testCase.WorkDir, glued{1}), '[redacted_url]');
        end
        function urlInAQueryLeavesWithTheQuery(testCase)
            testCase.verifyEqual(optiprofiler_internal.EvalReport.credentialFreeUrls( ...
                'https://example.org/x?next=https://u2:CANARYpw@example.org/y'), ...
                'https://example.org/x [url parameters removed]');
        end
        function longTextWithoutAUrlIsHandledInLinearTime(testCase)
            text = repmat('a', 1, 200000);
            start = tic;
            testCase.verifyEqual(numel(optiprofiler_internal.EvalReport.credentialFreeUrls(text)), 200000);
            testCase.verifyEqual(numel(TestReportTextSafety.encoded(text)), 256);
            testCase.verifyLessThan(toc(start), 5);
        end
        function absolutePathInsideATextIsReplaced(testCase)
            testCase.verifyEqual(TestReportTextSafety.encoded('built in /home/someone/x today'), ...
                'built in [redacted_absolute_path] today');
            testCase.verifyEqual(TestReportTextSafety.encoded('see C:\Users\someone\x'), 'see [redacted_absolute_path]');
        end
        function namesThatDifferOnlyInARedactedPartStayApart(testCase)
            first = optiprofiler_internal.EvalReport.problemId('user', 'ssh://user:one@example.org/a', 'primary');
            second = optiprofiler_internal.EvalReport.problemId('user', 'ssh://user:two@example.org/a', 'primary');
            testCase.verifyNotEqual(first, second);
            testCase.verifyFalse(contains(first, 'one'));
            testCase.verifyEqual(optiprofiler_internal.EvalReport.problemId('user', 'ROSENBR', 'primary'), ...
                jsonencode({'user', 'ROSENBR', 'primary'}));
        end
        function textWithoutAUrlIsNotTreatedAsOne(testCase)
            for text = {'not a url: a:b', 'mailto:someone@example.org'}
                testCase.verifyEqual(TestReportTextSafety.encoded(text{1}), text{1});
            end
            testCase.verifyEqual(TestReportTextSafety.encoded('C://Users/someone/x'), '[absolute_path_omitted]');
        end
        function sensitiveOptionNamesAreOmittedAtAnyDepth(testCase)
            for key = {'password', 'PassWord', 'user_passwd', 'ssh_passphrase', 'pwd', 'userPwd', 'client_secret', ...
                    'SECRET_KEY', 'aws_secret_access_key', 'credentials', 'Authorization', 'authentication', 'auth', ...
                    'http_auth', 'basicAuth', 'oauth', 'auth_value', 'bearer', 'bearer_value', 'api_key', 'apikey', ...
                    'API_KEY', 'private_key', 'privateKey', 'session_cookie', 'token', 'auth_token', 'authToken', ...
                    'access_token', 'sessiontoken', 'token_value'}
                value = TestReportTextSafety.encoded(struct(key{1}, 'do-not-copy', 'kept', 1));
                testCase.verifyEqual(value, struct('kept', 1), key{1});
                value = TestReportTextSafety.encoded(struct('outer', struct('inner', struct(key{1}, 'do-not-copy', 'kept', 1))));
                testCase.verifyEqual(value.outer.inner, struct('kept', 1), key{1});
            end
        end
        function ordinaryOptionNamesAreKept(testCase)
            for key = {'max_tokens', 'tokenizer', 'tokens', 'author', 'authors', 'authenticated', 'keyword', 'noise_level'}
                value = TestReportTextSafety.encoded(struct(key{1}, 3));
                testCase.verifyEqual(value, struct(key{1}, 3), key{1});
            end
            % A machine path is a privacy rule of its own, not a secret.
            value = TestReportTextSafety.encoded(struct('report_path', 'x', 'kept', 1));
            testCase.verifyEqual(value, struct('kept', 1));
        end
        function reportFilesHoldNoCredentialOfAnyUrlScheme(testCase)
            canary = 'CANARYpw';
            url = @(scheme) sprintf('%s://user:%s@example.org/repo.git', scheme, canary);
            reference = struct('fun', 0, 'maxcv', 0, 'kind', 'feasible', 'source', url('ssh'));
            options = struct('score_only', true, 'silent', true, 'n_jobs', 1, 'max_eval_factor', 2);
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, 'reference', reference, ...
                'name', ['problem from ', url('git+ssh')]));
            options.solver_names = {url('HTTPS'), url('git+https')};
            options.report_path = fullfile(testCase.WorkDir, 'report.json');
            benchmark({@TestReportTextSafety.finish, @TestReportTextSafety.stay}, options);
            files = dir(fullfile(testCase.WorkDir, '**', '*'));
            files = files(~[files.isdir]);
            testCase.verifyGreaterThanOrEqual(numel(files), 2);
            for k = 1:numel(files)
                text = fileread(fullfile(files(k).folder, files(k).name));
                testCase.verifyFalse(contains(text, canary), files(k).name);
            end
        end
        function savedReportAndArtifactNamesHoldNoCredential(testCase)
            canary = 'CANARYpw';
            options = struct('silent', true, 'n_jobs', 1, 'max_eval_factor', 2, 'draw_hist_plots', 'sequential');
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, ...
                'name', sprintf('problem from git+ssh://user:%s@example.org/repo.git', canary)));
            options.solver_names = {sprintf('ssh://u:%s@h', canary), 'b'};
            options.savepath = fullfile(testCase.WorkDir, 'out');
            options.report_path = fullfile(testCase.WorkDir, 'report.json');
            benchmark({@TestReportTextSafety.finish, @TestReportTextSafety.stay}, options);
            files = dir(fullfile(options.savepath, '**', '*'));
            testCase.verifyTrue(any(endsWith({files.name}, '.pdf')));
            testCase.verifyFalse(any(contains({files.name}, canary)));
            for name = {'report.json', 'report.plot_data.json'}
                testCase.verifyFalse(contains(fileread(fullfile(testCase.WorkDir, name{1})), canary), name{1});
            end
        end
        function failedRenderingScopeHoldsNoCredential(testCase)
            canary = 'CANARYpw';
            name = sprintf('ssh://user:%s@example.org/p', canary);
            target = fullfile(testCase.WorkDir, 'scope.json');
            report = optiprofiler_internal.EvalReport(target, struct(), @TestReportTextSafety.replaceFile);
            report.selection('user', {name}, 'primary');
            report.loadFailed('user', name, 'primary', 'Fixture:Load');
            report.recordRendering('failed', struct('library', 'user', 'name', name));
            report.finish();
            testCase.verifyFalse(contains(fileread(target), canary));
        end
        function longReferenceSourceIsMarkedWhenClipped(testCase)
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                ['best known value reported by ', repmat('y', 1, 300)]);
            testCase.verifyEqual(numel(reported), 256);
            testCase.verifyTrue(endsWith(reported, '...[truncated]'));
        end
    end
end
