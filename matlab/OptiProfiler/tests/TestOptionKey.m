classdef TestFeatureOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('FeatureOptionKey');
            expectedValues = {'distribution', 'modifier', 'n_runs', 'order', 'parameter', 'rate_error', 'rate_nan', 'significant_digits', 'type'};

            for i = 1:numel(enumValues)
                testCase.verifyEqual(enumValues(i).value, expectedValues{i});
            end
        end

    end

end