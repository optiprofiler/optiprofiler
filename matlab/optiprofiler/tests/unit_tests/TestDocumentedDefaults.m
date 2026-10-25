classdef TestDocumentedDefaults < matlab.unittest.TestCase
%TESTDOCUMENTEDDEFAULTS The benchmark help states the defaults the code applies.
    properties
        OldCwd
        OldPath
    end
    methods (TestMethodSetup)
        function setupState(testCase)
            testCase.OldCwd = pwd;
            testCase.OldPath = path;
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.OldCwd);
            path(testCase.OldPath);
        end
    end
    methods (Static)
        function x = stay(fun, x0, varargin)
            fun(x0);
            x = x0;
        end
    end
    methods (Test)
        function drawHistPlotsHelpMatchesDefault(testCase)
            text = help('benchmark');
            stated = regexp(text, 'draw_hist_plots:.*?Default is ''(\w+)''', 'tokens', 'once');
            testCase.assertNotEmpty(stated, 'The help text must state the draw_hist_plots default.');
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            options = getDefaultProfileOptions({@TestDocumentedDefaults.stay, @TestDocumentedDefaults.stay}, ...
                Feature('plain'), struct());
            testCase.verifyEqual(stated{1}, options.draw_hist_plots);
        end
    end
end
