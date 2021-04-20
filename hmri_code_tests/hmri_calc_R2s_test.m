% Unit tests implemented:
%   Multiple or single contrast weighting
%   1, 2 and 3D datasets tests, via permutation along dimensions
%   Zeros as input will return NaNs
%   Higher tolerance with noise added
%   Checks that input must be a structure

classdef hmri_calc_R2s_test < matlab.unittest.TestCase
    properties (TestParameter)
        sizes1 = {1,10,100};
        sizes2 = {1,10,100};
        sizes3 = {1,10,100};
        tolerance = {1e-9};
        noiseTol = {0.05};
    end
    properties
        rnd = [];
    end
    
    methods (Test)
        
        %% Test Functions
        function testSingleContrast(testCase,sizes1,sizes2,sizes3,tolerance)
            
            % Test on 3D simulated data that the calculated R2* is within a defined
            % tolerance given a single contrast input.
            % Single output (only) case is tested.
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            TEs=(2:2.5:20)*1e-3; % s
            signal_TE0=2000*rand(dims)+500;
            
            % Create signal decay given R2s and TEs:
            signal=hmri_calc_R2s_test.decaySignal(signal_TE0,TEs,R2s);
            
            % Check that it works without outputting the extrapolated values
            R2sEst=hmri_calc_R2s(struct('data',signal,'TE',TEs));
            
            % Check every element is less than pre-defined threshold, here 1e-9 %TODO
            % -rationale for threshold.
            assertEqual(testCase,R2sEst,R2s,'AbsTol',tolerance)
        end
        
       
        function testMultipleContrast(testCase,sizes1,sizes2,sizes3,tolerance)
            
            % Test on 3D simulated data that the calculated R2* and contrast-specific
            % intercepts are within a defined tolerance given multiple contrast input.
            % Two outputs case is tested.
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            TEs1=(2:2.5:20)*1e-3; % s
            signal1_TE0=2000*rand(dims)+500; % [500 2500]
            signal1=hmri_calc_R2s_test.decaySignal(signal1_TE0,TEs1,R2s);
            
            % First four TEs for second contrast
            TEs2=TEs1(1:4); % s
            signal2_TE0=1000*rand(dims)+100; % [100 1100]
            signal2=hmri_calc_R2s_test.decaySignal(signal2_TE0,TEs2,R2s);
            
            [R2sEst,extrapolated]=hmri_calc_R2s([struct('data',signal1,'TE',TEs1),struct('data',signal2,'TE',TEs2)]);
            
            assertEqual(testCase,R2sEst,R2s,'AbsTol',tolerance)
            assertEqual(testCase,extrapolated{1},signal1_TE0,'AbsTol',tolerance)
            assertEqual(testCase,extrapolated{2},signal2_TE0,'AbsTol',tolerance)
            
        end
        
        function testSingleContrast1D(testCase,sizes1,tolerance)
            
            % Test on 1D simulated data that the calculated R2* and intercept
            % are within a defined tolerance given single contrast input.
            
            dims=[sizes1,1];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            TEs=(2:2.5:20)*1e-3; % s
            signal_TE0=2000*rand(dims)+500;
            signal=hmri_calc_R2s_test.decaySignal(signal_TE0,TEs,R2s);
            
            [R2sEst,extrapolated]=hmri_calc_R2s(struct('data',signal,'TE',TEs));
            
            assertEqual(testCase,R2sEst,R2s,'AbsTol',tolerance)
            assertEqual(testCase,extrapolated{1},signal_TE0,'AbsTol',tolerance)
            
        end
        
        function testZero2DInputs(testCase,sizes1,sizes2,sizes3)
            
            % Test on 2D simulated data that the calculated R2* is NaN for a single
            % contrast input where the signal is actually zero.
            
            dims=[sizes1,sizes2,sizes3]; % last dim is echoes
            TEs=(2:2.5:20)*1e-3; % s
            signal=zeros([dims,length(TEs)]);
            
            % Check that it fails when zeros are input
            R2sEst=hmri_calc_R2s(struct('data',signal,'TE',TEs));
            
            % Check NaNs are returned
            assertTrue(testCase,all(isnan(R2sEst(:))))
            
        end
        
        function testNonStructInput(testCase)
            
            % Check that an error is thrown if the input isn't a struct
            assertError(testCase, @() hmri_calc_R2s(zeros([32,56,8])), 'hmri:structError');
            
        end
        
        function testNoise(testCase,sizes1,sizes2,sizes3,noiseTol)
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            % If max TE is too long, then estimation becomes inaccurate because of the
            % exponential decrease in SNR with TE
            TEs=(2:2.5:20)*1e-3; % s
            w_TE0=2000*rand(dims)+500;
            w=hmri_calc_R2s_test.decaySignal(w_TE0,TEs,R2s);
            
            % SNR = w_TE0
            wN=w+randn(size(w));
            
            [R2sEst,extrapolated]=hmri_calc_R2s(struct('data',wN,'TE',TEs));
            
            % Check relative error less than 5%
            assertEqual(testCase,R2sEst,R2s,'RelTol',noiseTol)
            assertEqual(testCase,extrapolated{1},w_TE0,'RelTol',noiseTol)
            
        end
        
        function testMultipleContrastNoise(testCase,sizes1,sizes2,sizes3,noiseTol)
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            % If max TE is too long, then estimation becomes inaccurate because of the
            % exponential decrease in SNR with TE
            TEs1=(2:2.5:20)*1e-3; % s
            
            w1_TE0=2000*rand(dims)+500;
            w1=hmri_calc_R2s_test.decaySignal(w1_TE0,TEs1,R2s);
            w1N=w1+randn(size(w1));
            
            TEs2=TEs1(1:6); % s
            w2_TE0=2000*rand(dims)+500;
            w2=hmri_calc_R2s_test.decaySignal(w2_TE0,TEs2,R2s);
            w2N=w2+randn(size(w2));
            
            [R2sEst,extrapolated]=hmri_calc_R2s([struct('data',w1N,'TE',TEs1),struct('data',w2N,'TE',TEs2)]);
            
            % Check relative error less than 5%
            assertEqual(testCase,R2sEst,R2s,'RelTol',noiseTol)
            assertEqual(testCase,extrapolated{1},w1_TE0,'RelTol',noiseTol)
            assertEqual(testCase,extrapolated{2},w2_TE0,'RelTol',noiseTol)
            
        end
    end
    
    %% Test Setup and Teardown Functions
    
    methods(TestMethodSetup)
        function seedRandomNumberGenerator(testCase) %#ok<MANU>
            % set randomness method and seed for reproducability
            %
            % it would be best to use a separate radom stream for each test
            % for simplicity we set the global stream
            % this will seed the RNG separately on every worker
            rng(319, 'twister');
        end
    end

    methods(TestClassSetup)
        function enableRngLegacyWarnings(testCase) %#ok<MANU>
            % make sure the use oflegacy code that breaks random number
            % handling is going to generate a warning
            warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
            warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState') 
        end
    end
    
    %% Auxilary Functions
    
    methods(Static)
        
        function w_TEs=decaySignal(w_TE0,TEs,R2s)
            
            dims=size(w_TE0);
            
            % Account for 1D case
            if (length(dims)==2)&&(dims(2)==1), dims=dims(1); end
            
            TEs=reshape(TEs,[ones(1,length(dims)),length(TEs)]);
            
            w_TEs=w_TE0.*exp(-R2s.*TEs);
            
        end

    end
end
