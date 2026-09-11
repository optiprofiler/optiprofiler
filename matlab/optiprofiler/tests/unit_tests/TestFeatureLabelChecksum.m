classdef TestFeatureLabelChecksum < matlab.unittest.TestCase
% The display-label checksum is standard CRC-32, never a scientific seed.
    properties (Access = private)
        GetDefaults
    end
    methods (TestMethodSetup)
        function useRealDefaultBuilder(testCase)
            original = pwd;
            cleanup = onCleanup(@() cd(original)); %#ok<NASGU>
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            testCase.GetDefaults = @getDefaultProfileOptions;
        end
    end
    methods (Test)
        function standardUtf8Vectors(testCase)
            checksum = @optiprofiler_internal.featureLabelChecksum;
            testCase.verifyEqual(checksum(''), '00000000');
            testCase.verifyEqual(checksum('123456789'), 'cbf43926');
            % Independent Python zlib.crc32 UTF-8 reference, computed on syu.
            testCase.verifyEqual(checksum(['OptiProfiler ', char(955)]), 'b35ea718');
        end

        function generatedBoundaryAndFullLabelAreRetained(testCase)
            names = {{'custom', 'custom', 'custom', 'truncated', 'noisy'}, ...
                {'custom', 'custom', 'custom', 'custom', 'random_nan', 'random_nan'}, ...
                {'custom', 'custom', 'custom', 'truncated', 'truncated', 'random_nan'}};
            full = {'custom__custom__custom__truncated_6__noisy_0.001_mixed_gaussian', ...
                'custom__custom__custom__custom__random_nan_0.05__random_nan_0.05', ...
                'custom__custom__custom__truncated_6__truncated_6__random_nan_0.05'};
            get_defaults = testCase.GetDefaults;
            options = struct('n_jobs', 1, 'silent', true, 'score_only', true);
            for i = 1:3
                [effective, retained] = get_defaults({}, Feature(names{i}), options);
                testCase.verifyEqual(length(full{i}), 62 + i);
                testCase.verifyEqual(retained, full{i});
                if i < 3
                    testCase.verifyEqual(effective.feature_stamp, full{i});
                else
                    testCase.verifyEqual(effective.feature_stamp, ...
                        [full{i}(1:55), '_', optiprofiler_internal.featureLabelChecksum(full{i})]);
                    testCase.verifyEqual(length(effective.feature_stamp), 64);
                end
            end
            reversed = fliplr(names{3});
            [reverse_options, reverse_full] = get_defaults({}, Feature(reversed), options);
            testCase.verifyNotEqual(reverse_full, full{3});
            testCase.verifyNotEqual(reverse_options.feature_stamp, effective.feature_stamp);
        end

        function atomicAndExplicitLabelsStayUnchanged(testCase)
            get_defaults = testCase.GetDefaults;
            options = struct('n_jobs', 1, 'silent', true, 'score_only', true);
            [atomic, full] = get_defaults({}, Feature('noisy'), options);
            testCase.verifyEqual(atomic.feature_stamp, 'noisy_0.001_mixed_gaussian');
            testCase.verifyEqual(full, atomic.feature_stamp);
            % This is a literal user override, not a generated declaration.
            options.feature_stamp = repmat('x', 1, 80);
            [explicit, full] = get_defaults({}, Feature({'noisy', 'noisy'}), options);
            testCase.verifyEqual(explicit.feature_stamp, options.feature_stamp);
            testCase.verifyEqual(full, options.feature_stamp);
        end
    end
end
