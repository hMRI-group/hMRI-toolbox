% Unit tests based on input/output:
%   hmri_coreg
%   

classdef hmri_coreg_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.

    end
    properties
        % Parameters independent of TestParameter
        ref_file
        ref
        src_file
        src
        coregTol = 0.005
    end
    
    methods (Test)
        
        %% Test Functions
        function test_hmri_coreg(testCase)
            
            % Test that coregistration works on sample data.
            % Low resolution data used to speed up the test.
            
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
            assertEqual(testCase,x_est,x_app(1:length(x_est)),'RelTol',testCase.coregTol)
            
        end
            
 
    end
    
    %% Test Setup and Teardown Functions
    
    methods(TestMethodSetup)
        % These methods are run before each test
        
        function getData(testCase)
            % These methods are run when instantiating the class
            
            % Use low resolution brain image
            ut_data_dir = [fileparts(which('hmri_test_utils')) filesep 'example_data'];
            field_map_1 = [ut_data_dir filesep 'field_map_1.nii'];
            assert(logical(exist(field_map_1,'file')),'Could not find\n\t%s.\nPlease run "hmri_get_ut_data" to download the data',field_map_1)
            
            % Use temporary directory which is deleted after tests have run
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tempFixture = testCase.applyFixture(TemporaryFolderFixture);
            
            % Copy brain image to temporary directory; source and reference
            % images (for coregistration) are initially identical
            testCase.ref_file = [tempFixture.Folder filesep 'field_map_1.nii'];
            testCase.src_file = [tempFixture.Folder filesep 'field_map_1_copy.nii'];
            
            copyfile(field_map_1, testCase.ref_file);
            copyfile(field_map_1, testCase.src_file);
            
            testCase.ref = spm_vol(testCase.ref_file);
            testCase.src = spm_vol(testCase.src_file);
        end
        
    end

    methods(TestClassSetup)
        % These methods are run before each test suite (i.e. before all 
        % tests of this class are run)
  
    end
    
    methods(TestMethodTeardown)
        
    end
    
    methods(TestClassTeardown)

    end
 
end
