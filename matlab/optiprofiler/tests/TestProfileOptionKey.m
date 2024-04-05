classdef TestProfileOptionKey < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('ProfileOptionKey');
            expectedValues = {'N_JOBS'};

            diff1 = setdiff(enumValues, expectedValues);
            diff2 = setdiff(expectedValues, enumValues);
            testCase.verifyEmpty(diff1, "The enumeration values are not correctly assigned");
            testCase.verifyEmpty(diff2, "The enumeration values are not correctly assigned");
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected
            testCase.verifyEqual(ProfileOptionKey.N_JOBS.value, 'n_jobs');
        end
        
    end

end