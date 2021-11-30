% Unit tests implemented:
%

classdef hmri_read_vols_test < matlab.unittest.TestCase
    properties (TestParameter)
        V = {struct}
        VG = {struct};
        interp = {1};
        p = {1};
        x = {1}
    end
    
    methods (Test)
        
        %% Test Functions
        function testNonStructInput(testCase, V, VG, p, interp, x)
            
            % Check that an error is thrown if the input isn't a struct
            assertError(testCase, @() hmri_read_vols('myDummyUTNonStruct.nii', VG, p, interp, x), 'hmri:structError');
            assertError(testCase, @() hmri_read_vols(V, 'myDummyUTNonStruct.nii', p, interp, x), 'hmri:structError');
            assertError(testCase, @() hmri_read_vols(V, VG, p, 130, x), 'hmri:inputError');
            assertError(testCase, @() hmri_read_vols(V, VG, 'test', interp, x), 'hmri:typeError');
            assertError(testCase, @() hmri_read_vols(V, VG, p, interp, ones(8,7)), 'hmri:typeError');
            
        end
    end
end