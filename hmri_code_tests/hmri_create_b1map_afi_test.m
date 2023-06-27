% Unit tests for AFI B1 map creation

classdef hmri_create_b1map_afi_test < matlab.unittest.TestCase
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
        function test_afi(testCase)
            
            % Test that AFI B1 mapping works with synthetic data like from 
            % BIDS and hMRI toolbox DICOM conversion

            V=testCase.phantomStruct;

            %% Create synthetic AFI data from phantom B1 map
            B1=spm_read_vols(V)/100;
            
            % AFI acquisition parameters: different from defaults
            alphanom=55;            
            TR2TR1ratio=10;

            TR1=30; % ms
            TR2=TR1*TR2TR1ratio;

            % Tissue parameters
            A=1e3;
            R1=1/1e3; % 1/ms
            
            % TR1
            S1=abs(A*hmri_test_utils.dualTRernstd(B1*alphanom,TR1,TR2,R1)); 
            fname1=spm_file(testCase.phantomStruct.fname,'suffix','_AFI_TR1'); 
            V.fname=fname1;
            V.descrip=sprintf('Synthetic AFI data for hMRI toolbox testing TR=%dms/TE=%dms/FA=%ddeg',TR1,0,alphanom);
            spm_write_vol(V,S1);

            % TR2
            S2=abs(A*hmri_test_utils.dualTRernstd(B1*alphanom,TR2,TR1,R1));
            fname2=spm_file(testCase.phantomStruct.fname,'suffix','_AFI_TR2'); 
            V.fname=fname2;
            V.descrip=sprintf('Synthetic AFI data for hMRI toolbox testing TR=%dms/TE=%dms/FA=%ddeg',TR2,0,alphanom);
            spm_write_vol(V,S2);

            %% Create SPM batch job to run
            config=fullfile(fileparts(mfilename('fullpath')),'testconfig','hmri_b1_local_defaults_test_afi_bids.m');

            inputs{1} = {fname1; fname2}; % Create B1 map: B1 input - cfg_files
            inputs{2} = {config};         % Create B1 map: Customised B1 defaults file - cfg_files

            matlabbatch{1}.spm.tools.hmri.create_B1.subj.output.indir = 'yes';
            matlabbatch{1}.spm.tools.hmri.create_B1.subj.b1_type.i3D_AFI.b1input = '<UNDEFINED>';
            matlabbatch{1}.spm.tools.hmri.create_B1.subj.b1_type.i3D_AFI.b1parameters.b1defaults = '<UNDEFINED>';
            matlabbatch{1}.spm.tools.hmri.create_B1.subj.popup = false;
            
            %% Run SPM batch job and evaluate result with bids sidecar file
            bids.FlipAngle=alphanom;

            bids.RepetitionTimeExcitation=TR1;
            f=fopen(spm_file(fname1,'ext','.json'),'w');
            fprintf(f,jsonencode(bids));
            fclose(f);

            bids.RepetitionTimeExcitation=TR2;
            f=fopen(spm_file(fname2,'ext','.json'),'w');
            fprintf(f,jsonencode(bids));
            fclose(f);

            spm('defaults', 'FMRI');
            out = spm_jobman('run', matlabbatch, inputs{:});
            
            B1est = spm_read_vols(spm_vol(out{1}.subj.B1map{1}))/100;
            assertEqual(testCase, B1est, B1, 'RelTol', 1e-2)

            %% Run SPM batch job and evaluate result with hmri-like sidecar file
            hmri.acqpar.FlipAngle=alphanom;
            hmri.acqpar.RepetitionTime=TR1; % same TR for both stored in DICOM
            hmri.acqpar.alTR=[TR1,TR2];

            f=fopen(spm_file(fname1,'ext','.json'),'w');
            fprintf(f,jsonencode(hmri));
            fclose(f);

            f=fopen(spm_file(fname2,'ext','.json'),'w');
            fprintf(f,jsonencode(hmri));
            fclose(f);

            spm('defaults', 'FMRI');
            out = spm_jobman('run', matlabbatch, inputs{:});
            
            B1est = spm_read_vols(spm_vol(out{1}.subj.B1map{1}))/100;
            assertEqual(testCase, B1est, B1, 'RelTol', 1e-2)

            %% Test that appropriate warning is given if images given in wrong order
            inputs_wrong = inputs;
            inputs_wrong{1} = {fname2; fname1}; % Create B1 map: B1 input - cfg_files

            spm('defaults', 'FMRI');
            verifyWarning(testCase,@() spm_jobman('run', matlabbatch, inputs_wrong{:}),'hmri:afiTooManyImag')
        end
        
        
    end
    
    %% Test Setup and Teardown Functions
    
    methods(TestMethodSetup)
        % These methods are run before each test
        
    end
    
    methods(TestClassSetup)
        % These methods are run when instantiating the class
        function resetToDefaults(testCase)
            hmri_b1_standard_defaults;
            testCase.addTeardown(@hmri_b1_standard_defaults)
        end
        
        % Create synthetic phantom
        function generatePhantom(testCase)
            % Temporary fixtures are deleted after tests complete
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tempFixture = testCase.applyFixture(TemporaryFolderFixture);
            
            [testCase.phantomData,testCase.phantomStruct]=...
                hmri_test_utils.makePhantom_B1([30,30,30],[70,120],fullfile(tempFixture.Folder,'phantomB1.nii'));
        end
    end
    
end
