classdef TestNativeGraphicsTempPath < matlab.unittest.TestCase
%TESTNATIVEGRAPHICSTEMPPATH Verified Linux renderer boundary, without globals.
    properties (Access=private)
        CheckPath
        CheckGraphics
    end
    methods (TestMethodSetup)
        function findPrivatePolicy(testCase)
            old=pwd; cleanup=onCleanup(@() cd(old)); %#ok<NASGU>
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            testCase.CheckPath=@nativeGraphicsTempPathIsUnsafe;
            testCase.CheckGraphics=@hasNativeGraphics;
        end
    end
    methods (Test)
        function actualEnvironmentReachesGraphicsPolicy(testCase)
            testCase.assumeEqual(computer('arch'), 'glnxa64');
            testCase.assumeEqual(version('-release'), '2026a');
            % setenv(name, getenv(name)) would turn an unset TMPDIR into an
            % explicitly empty one, which disables native export for the rest
            % of the session.
            was_set = isenv('TMPDIR');
            original = getenv('TMPDIR');
            testCase.addTeardown(@() restoreTmpdir(was_set, original));
            unsafe = ['/tmp/', repmat('a', 1, 57)];
            setenv('TMPDIR', unsafe);
            testCase.verifyFalse(testCase.CheckGraphics());
            testCase.verifyEqual(getenv('TMPDIR'), unsafe);
        end
        function restoringAnUnsetTmpdirLeavesItUnset(testCase)
            testCase.assumeTrue(exist('isenv') > 0 && exist('unsetenv') > 0, 'isenv/unsetenv need R2022b or later.'); %#ok<EXIST>
            was_set = isenv('TMPDIR');
            original = getenv('TMPDIR');
            testCase.addTeardown(@() restoreTmpdir(was_set, original));
            unsetenv('TMPDIR');
            setenv('TMPDIR', '/tmp/x');
            restoreTmpdir(false, '');
            testCase.verifyFalse(isenv('TMPDIR'));
            setenv('TMPDIR', '/tmp/x');
            restoreTmpdir(true, '/tmp/y');
            testCase.verifyEqual(getenv('TMPDIR'), '/tmp/y');
        end
        function linuxR2026aByteBoundary(testCase)
            testCase.verifyFalse(testCase.CheckPath(['/tmp/',repmat('a',1,56)],'glnxa64','2026a'));
            testCase.verifyTrue(testCase.CheckPath(['/tmp/',repmat('a',1,57)],'glnxa64','2026a'));
        end
        function countsUtf8BytesNotCharacters(testCase)
            % U+00E9 is two UTF-8 bytes but one MATLAB character.
            safe=['/tmp/',repmat('a',1,54),char(233)];
            unsafe=['/tmp/',repmat('a',1,55),char(233)];
            testCase.verifyEqual(numel(unicode2native(safe,'UTF-8')),61);
            testCase.verifyEqual(numel(unicode2native(unsafe,'UTF-8')),62);
            testCase.verifyFalse(testCase.CheckPath(safe,'glnxa64','2026a'));
            testCase.verifyTrue(testCase.CheckPath(unsafe,'glnxa64','2026a'));
        end
        function trailingSeparatorsDoNotChangeBoundary(testCase)
            testCase.verifyFalse(testCase.CheckPath(['/tmp/',repmat('a',1,56),'///'],'glnxa64','2026a'));
            testCase.verifyTrue(testCase.CheckPath(['/tmp/',repmat('a',1,57),'///'],'glnxa64','2026a'));
            testCase.verifyFalse(testCase.CheckPath('////','glnxa64','2026a'));
        end
        function relativePathsUseConservativePolicy(testCase)
            for value={'relative','./tmp','../tmp','~/.tmp'}
                [unsafe,reason]=testCase.CheckPath(value{1},'glnxa64','2026a');
                testCase.verifyTrue(unsafe);
                testCase.verifySubstring(reason,'relative');
            end
        end
        function shortAbsolutePathsRemainNative(testCase)
            for value={'/','/tmp','/tmp/with space'}
                [unsafe,reason]=testCase.CheckPath(value{1},'glnxa64','2026a');
                testCase.verifyFalse(unsafe);
                testCase.verifyEmpty(reason);
            end
        end
        function unsetAndExplicitlyEmptyAreDifferent(testCase)
            [unsafe,reason]=testCase.CheckPath('','glnxa64','2026a',false);
            testCase.verifyFalse(unsafe); testCase.verifyEmpty(reason);
            [unsafe,reason]=testCase.CheckPath('','glnxa64','2026a',true);
            testCase.verifyTrue(unsafe); testCase.verifySubstring(reason,'explicitly set');
        end
        function unknownEmptyPresenceUsesConservativePolicy(testCase)
            [unsafe,reason]=testCase.CheckPath('','glnxa64','2026a',[]);
            testCase.verifyTrue(unsafe);
            testCase.verifySubstring(reason,'could not be established');
        end
        function otherReleasesRemainUnchanged(testCase)
            for release={'2023b','2025b','2026b'}
                testCase.verifyFalse(testCase.CheckPath(['/',repmat('x',1,200)],'glnxa64',release{1}));
                testCase.verifyFalse(testCase.CheckPath('relative','glnxa64',release{1}));
                testCase.verifyFalse(testCase.CheckPath('','glnxa64',release{1},true));
            end
        end
        function otherPlatformsRemainUnchanged(testCase)
            for platform={'win64','maci64','maca64'}
                testCase.verifyFalse(testCase.CheckPath(['/',repmat('x',1,200)],platform{1},'2026a'));
                testCase.verifyFalse(testCase.CheckPath('relative',platform{1},'2026a'));
                testCase.verifyFalse(testCase.CheckPath('',platform{1},'2026a',true));
            end
        end
    end
end

function restoreTmpdir(was_set, value)
    if was_set
        setenv('TMPDIR', value);
    else
        unsetenv('TMPDIR');
    end
end
