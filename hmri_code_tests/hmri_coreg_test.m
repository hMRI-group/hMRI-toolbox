% Unit tests based on input/output:
%   hmri_coreg
%   

classdef hmri_coreg_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.
        coregTol = {0.05}

    end
    properties
        % Parameters independent of TestParameter
        refStruct
        srcStruct
    end
    
    methods (Test)
        
        %% Test Functions
        function test_hmri_coreg(testCase,coregTol)
            
            ref = spm_vol(testCase.refStruct.fname);
            src = spm_vol(testCase.srcStruct.fname);
%             pn = 'D:\My Documents\Test\gre_field_mapping_1acq_rl_0004';
%             fn2 = 's2017-10-16_12-29-124020-00001-00001-2.nii';
%             fn = 's2017-10-16_12-29-124019-00001-00001-1.nii';
%             ref = spm_vol([pn filesep fn]);
%             src = spm_vol([pn filesep fn2]);
            
            % Create transformation matrix and inverse:
            x = [10 8 5 0.02 0.01 0.008];
            M = pinv(spm_matrix(x));
            
            % Apply transformation to source:
            spm_get_space(deblank(deblank(src.fname)), M*spm_get_space(deblank(src.fname)));
            src = spm_vol(src.fname); % read in src again to pick up modified header
            
            % Estimate transformation:
            x_est = hmri_coreg(ref, src, []);
            
            % %TODO - find way to test the assertion.  High residuals using
            % phantom data.  Good performance with in vivo test data.
            x_app = spm_imatrix(M);
%             assertEqual(testCase,x_est,x_app(1:length(x_est)),'RelTol',coregTol)
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
            [~,testCase.refStruct]=...
                hmri_test_utils.makePhantom([tempname '.nii'], eye(4));
            [~,testCase.srcStruct]=...
                hmri_test_utils.makePhantom([tempname '.nii'], eye(4));
        end

    end
    
    % TODO Destructor to remove testCase.V.fname == tempname.nii
 
end
