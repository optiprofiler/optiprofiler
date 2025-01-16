classdef TestNoiseType < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            
            enumValues = enumeration('NoiseType');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'ABSOLUTE'; 'RELATIVE'; 'MIXED'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected
            
            clear obj;
            clear NoiseType;
            testCase.verifyEqual(NoiseType.ABSOLUTE.value, 'absolute');
            testCase.verifyEqual(NoiseType.RELATIVE.value, 'relative');
            testCase.verifyEqual(NoiseType.MIXED.value, 'mixed');
        end
        
    end

end