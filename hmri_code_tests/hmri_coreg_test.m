% Unit tests based on input/output:
%   hmri_coreg
%   

classdef hmri_coreg_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.
        coregTol = {0.001}

    end
    properties
        % Parameters independent of TestParameter
        src_file
        src
        ref
    end
    
    methods (Test)
        
        %% Test Functions
        function test_hmri_coreg(testCase,coregTol)
            
            % Create transformation matrix and inverse:
            x = [10 8 5 0.02 0.01 0.008];
            M = pinv(spm_matrix(x));
            
            % Apply transformation to source:
            spm_get_space(deblank(deblank(testCase.src.fname)), M*spm_get_space(deblank(testCase.src.fname)));
            testCase.src = spm_vol(testCase.src.fname); % read in src again to pick up modified header
            
            % Estimate transformation:
            x_est = hmri_coreg(testCase.ref, testCase.src, []);
            
            % Test equivalence after coregistration:
            x_app = spm_imatrix(M);
            assertEqual(testCase,x_est,x_app(1:length(x_est)),'RelTol',coregTol)
            
        end
            
 
    end
    
    %% Test Setup and Teardown Functions
    
    methods(TestMethodSetup)
        % These methods are run before each test
        
        function getData(testCase)
            % These methods are run when instantiating the class
            ut_data_dir = [fileparts(which('hmri_test_utils')) filesep 'example_data'];
            ref_file = [ut_data_dir filesep 'field_map_1.nii'];
            testCase.src_file = [ut_data_dir filesep 'field_map_1_copy.nii'];
            testCase.ref = spm_vol(ref_file);
            copyfile(ref_file, testCase.src_file);
            testCase.src = spm_vol(testCase.src_file);
        end
        
    end

    methods(TestClassSetup)
        % These methods are run before each test suite (i.e. before all 
        % tests of this class are run)
  
    end
    
    methods(TestMethodTeardown)
        
        function deleteTempData(testCase)
            % Destructor 
            delete(testCase.src_file)
            
            % Should also delete ps coreg output created by spm
        end
    end
    
    methods(TestClassTeardown)

    end
 
end
