classdef TestNativeGraphicsTempPath < matlab.unittest.TestCase
%TESTNATIVEGRAPHICSTEMPPATH Verified Linux renderer boundary, without globals.
    properties (Access=private)
        CheckPath
    end
    methods (TestMethodSetup)
        function findPrivatePolicy(testCase)
            old=pwd; cleanup=onCleanup(@() cd(old)); %#ok<NASGU>
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            testCase.CheckPath=@nativeGraphicsTempPathIsUnsafe;
        end
    end
    methods (Test)
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
