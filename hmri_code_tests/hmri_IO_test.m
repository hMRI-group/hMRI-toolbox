% Unit tests based on input/output:
%   hmri_create_nifti - to do
%   hmri_read_vols - to do
%   

classdef hmri_IO_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.

    end
    properties
        % Parameters independent of TestParameter
        phantomData
        phantomStruct
    end
    
    methods (Test)
        
        %% Test Functions
        function test_hmri_create_nifti(testCase)
            
            dt = [spm_type('int32'),spm_platform('bigend')]; % for nifti output
            Ni = hmri_create_nifti([tempname '.nii'], testCase.phantomStruct, dt, 'test_hmri_create_nifti');
            Ni.dat(:,:,:) = testCase.phantomData;
            
            data = spm_read_vols(spm_vol(Ni.dat.fname));
            refData = spm_read_vols(spm_vol(testCase.phantomStruct.fname));
            
            % Check every element of R2sEst is equal to R2s to within a 
            % pre-defined threshold, defined as object property. 
            % %TODO -rationale for threshold.
            assertEqual(testCase,data,refData)
        end
            
 
    end
    
    %% Test Setup and Teardown Functions
    
    methods(TestMethodSetup)
        % These methods are run before each test
        
    end

    methods(TestClassSetup)
        % These methods are run when instantiating the class
        
        % Create synthetic phantom
        function generatePhantom(testCase)
            [testCase.phantomData,testCase.phantomStruct]=...
                hmri_test_utils.makePhantom([tempname '.nii'], eye(4));
        end

    end
    
    % TODO Destructor to remove testCase.V.fname == tempname.nii
 
end
