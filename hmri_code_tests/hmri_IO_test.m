% Unit tests based on input/output

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
            
            % Test that hmri_create_nifti writes out the same data it is given
            
            dt = [spm_type('int32'),spm_platform('bigend')]; % for nifti output
            Ni = hmri_create_nifti([tempname '.nii'], testCase.phantomStruct, dt, 'test_hmri_create_nifti');
            Ni.dat(:,:,:) = testCase.phantomData;
            
            data = spm_read_vols(spm_vol(Ni.dat.fname));
            refData = spm_read_vols(spm_vol(testCase.phantomStruct.fname));
            
            assertEqual(testCase,data,refData)
        end
        
        function test_hmri_read_vols(testCase)
            
            % Test that hmri_read_vols correctly reads in data when no 
            % transformation of the header is needed
            
            V=testCase.phantomStruct;
            data=zeros(V.dim);
            for p=1:V.dim(end)
                data(:,:,p) = hmri_read_vols(V,V,p,0);
            end
            refData = double(testCase.phantomData);
            
            assertEqual(testCase,data,refData)
        end
        
        function test_hmri_read_vols_hdr_shift(testCase)
            
            % Test that hmri_read_vols reads in data correctly after the
            % transformation in the header has been altered by a known
            % amount
            
            VG=testCase.phantomStruct;
            
            % Create transformation matrix and inverse
            x = [10 8 5 0.02 0.01 0.008];
            M = pinv(spm_matrix(x));
            
            % Apply transformation to source header
            srcfname=fullfile(spm_file(VG.fname,'path'),'source.nii');
            copyfile(VG.fname,srcfname);
            spm_get_space(srcfname, M*spm_get_space(srcfname));
            V = spm_vol(srcfname); % read in src to pick up modified header
            
            % Undo the header transformation when reading in data
            data=zeros(VG.dim);
            for p=1:VG.dim(end)
                data(:,:,p) = hmri_read_vols(V,VG,p,0,spm_imatrix(M));
            end
            refData = double(testCase.phantomData);
            
            assertEqual(testCase,data,refData)
        end
        
        function test_hmri_read_vols_NonStructInput(testCase)
            
            % Test that hmri_read_vols fails in the expected way when it is
            % given the wrong type of input
            
            % Minimal inputs just for testing
            V = struct;
            VG = struct;
            interp = 1;
            p = 1;
            x = 1;
            
            % Check that an error is thrown if the input isn't a struct
            assertError(testCase, @() hmri_read_vols('myDummyUTNonStruct.nii', VG, p, interp, x), 'hmri:structError');
            assertError(testCase, @() hmri_read_vols(V, 'myDummyUTNonStruct.nii', p, interp, x), 'hmri:structError');
            assertError(testCase, @() hmri_read_vols(V, VG, p, 130, x), 'hmri:inputError');
            assertError(testCase, @() hmri_read_vols(V, VG, 'test', interp, x), 'hmri:typeError');
            assertError(testCase, @() hmri_read_vols(V, VG, p, interp, ones(8,7)), 'hmri:typeError');
            
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
            % Temporary fixtures are deleted after tests complete
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tempFixture = testCase.applyFixture(TemporaryFolderFixture);
            
            [testCase.phantomData,testCase.phantomStruct]=...
                hmri_test_utils.makePhantom(fullfile(tempFixture.Folder,'phantom.nii'), eye(4));
        end
    end
    
end
