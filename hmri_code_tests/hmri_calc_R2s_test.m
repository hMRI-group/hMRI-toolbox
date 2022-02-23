% Unit tests implemented:
%   Multiple or single contrast weighting
%   1, 2 and 3D datasets tests, via permutation along dimensions
%   Zeros as input will return NaNs
%   Higher tolerance with noise added
%   Checks that input must be a structure

classdef hmri_calc_R2s_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.
        sizes1 = {1,10,100};
        sizes2 = {1,10};
        sizes3 = {1,10};
        fitmethod = {'OLS','WLS1','WLS3','NLLS_OLS'};
    end
    properties
        % Parameters independent of TestParameter
        tolerance = 1e-9; % absolute tolerance for numerical precision computations
        noiseTol = 0.05; % relative tolerance for testing the influence of noise
    end
    
    methods (Test)
        
        %% Test Functions
        function testSingleContrast(testCase,sizes1,sizes2,sizes3)
            
            % Test on 3D simulated data that the calculated R2* is within a defined
            % tolerance given a single contrast input.
            % Single output (only) case is tested.
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            TEs=(2:2.5:20)*1e-3; % s
            signal_TE0=2000*rand(dims)+500;
            
            % Create signal decay given R2s and TEs:
            signal=hmri_test_utils.createDecaySignal(signal_TE0,TEs,R2s);
            
            % Check that it works without outputting the extrapolated values
            R2sEst=hmri_calc_R2s(struct('data',signal,'TE',TEs),'OLS');
            
            % Check every element of R2sEst is equal to R2s to within a 
            % pre-defined threshold, defined as object property. 
            % %TODO -rationale for threshold.
            assertEqual(testCase,R2sEst,R2s,'AbsTol',testCase.tolerance)
        end
        
        function testMultipleContrast(testCase,sizes1,sizes2,sizes3)
            
            % Test on 3D simulated data that the calculated R2* and contrast-specific
            % intercepts are within a defined tolerance given multiple contrast input.
            % Two outputs case is tested.
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            TEs1=(2:2.5:20)*1e-3; % s
            signal1_TE0=2000*rand(dims)+500; % [500 2500]
            signal1=hmri_test_utils.createDecaySignal(signal1_TE0,TEs1,R2s);
            
            % First four TEs for second contrast
            TEs2=TEs1(1:4); % s
            signal2_TE0=1000*rand(dims)+100; % [100 1100]
            signal2=hmri_test_utils.createDecaySignal(signal2_TE0,TEs2,R2s);
            
            [R2sEst,extrapolated]=hmri_calc_R2s([struct('data',signal1,'TE',TEs1),struct('data',signal2,'TE',TEs2)],'OLS');
            
            assertEqual(testCase,R2sEst,R2s,'AbsTol',testCase.tolerance)
            assertEqual(testCase,extrapolated{1},signal1_TE0,'AbsTol',testCase.tolerance)
            assertEqual(testCase,extrapolated{2},signal2_TE0,'AbsTol',testCase.tolerance)
            
        end
        
        function testSingleContrast1D(testCase,sizes1)
            
            % Test on 1D simulated data that the calculated R2* and intercept
            % are within a defined tolerance given single contrast input.
            
            dims=[sizes1,1];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            TEs=(2:2.5:20)*1e-3; % s
            signal_TE0=2000*rand(dims)+500;
            signal=hmri_test_utils.createDecaySignal(signal_TE0,TEs,R2s);
            
            [R2sEst,extrapolated]=hmri_calc_R2s(struct('data',signal,'TE',TEs),'OLS');
            
            assertEqual(testCase,R2sEst,R2s,'AbsTol',testCase.tolerance)
            assertEqual(testCase,extrapolated{1},signal_TE0,'AbsTol',testCase.tolerance)
            
        end
        
        function testZero2DInputs(testCase,sizes1,sizes2,sizes3)
            
            % Test on 2D simulated data that the calculated R2* is NaN for a single
            % contrast input where the signal is actually zero.
            
            dims=[sizes1,sizes2,sizes3]; % last dim is echoes
            TEs=(2:2.5:20)*1e-3; % s
            signal=zeros([dims,length(TEs)]);
            
            % Check that it fails when zeros are input
            warning('off','hmri:zerosInInput')
            R2sEst=hmri_calc_R2s(struct('data',signal,'TE',TEs),'OLS');
            warning('on','hmri:zerosInInput')
            
            % Check NaNs are returned
            assertTrue(testCase,all(isnan(R2sEst(:))))
            
        end
        
        function testNonStructInput(testCase)
            
            % Check that an error is thrown if the input isn't a struct
            assertError(testCase, @() hmri_calc_R2s(zeros([32,56,8]),'OLS'), 'hmri:structError');
            
        end
        
        function testNoise(testCase,sizes1,sizes2,sizes3,fitmethod)
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            % If max TE is too long, then estimation becomes inaccurate because of the
            % exponential decrease in SNR with TE
            TEs=(2:2.5:20)*1e-3; % s
            w_TE0=2000*rand(dims)+500;
            w=hmri_test_utils.createDecaySignal(w_TE0,TEs,R2s);
            
            % SNR = w_TE0
            wN=abs(w+randn(size(w)));
            
            [R2sEst,extrapolated]=hmri_calc_R2s(struct('data',wN,'TE',TEs),fitmethod);
            
            % Check R2sEst and R2s are equal to within a relative error,
            % specified as a class property.
            assertEqual(testCase,R2sEst,R2s,'RelTol',testCase.noiseTol)
            assertEqual(testCase,extrapolated{1},w_TE0,'RelTol',testCase.noiseTol)
            
        end
        
        function testMultipleContrastNoise(testCase,sizes1,sizes2,sizes3,fitmethod)
            
            dims=[sizes1,sizes2,sizes3];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            % If max TE is too long, then estimation becomes inaccurate because of the
            % exponential decrease in SNR with TE
            TEs1=(2:2.5:20)*1e-3; % s
            
            w1_TE0=2000*rand(dims)+500;
            w1=hmri_test_utils.createDecaySignal(w1_TE0,TEs1,R2s);
            w1N=abs(w1+randn(size(w1)));
            
            TEs2=TEs1(1:6); % s
            w2_TE0=2000*rand(dims)+500;
            w2=hmri_test_utils.createDecaySignal(w2_TE0,TEs2,R2s);
            w2N=abs(w2+randn(size(w2)));
            
            [R2sEst,extrapolated]=hmri_calc_R2s([struct('data',w1N,'TE',TEs1),struct('data',w2N,'TE',TEs2)],fitmethod);
            
            % Check R2sEst and R2s are equal to within a relative error,
            % specified as a class property.
            assertEqual(testCase,R2sEst,R2s,'RelTol',testCase.noiseTol)
            assertEqual(testCase,extrapolated{1},w1_TE0,'RelTol',testCase.noiseTol)
            assertEqual(testCase,extrapolated{2},w2_TE0,'RelTol',testCase.noiseTol)
            
        end
        
        function testOLSvsWLS(testCase)
            
            dims=[100,10,10];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            % If max TE is too long, then estimation becomes inaccurate because of the
            % exponential decrease in SNR with TE
            TEs=(2:2.5:20)*1e-3; % s
            w_TE0=2000*rand(dims)+500;
            w=hmri_test_utils.createDecaySignal(w_TE0,TEs,R2s);
            
            % Add noise
            sigma=1; % sigma = 1 => SNR = w_TE0
            wN=abs(w+sigma*randn(size(w)));
            
            % Fit using OLS and WLS
            [R2sEstOLS,w_TE0_EstOLS]=hmri_calc_R2s(struct('data',wN,'TE',TEs),'OLS');
            [R2sEstWLS,w_TE0_EstWLS]=hmri_calc_R2s(struct('data',wN,'TE',TEs),'WLS1');
            
            % Compute theoretical covariance matrix of parameter estimates
            % https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation
            D=[TEs(:),ones(length(TEs),1)];
            errOLS=zeros([2,2,dims]);
            errWLS=zeros([2,2,dims]);
            for m=1:prod(dims)
                % Covariance matrix of additive noise after log transform
                M=diag((sigma./hmri_test_utils.createDecaySignal(w_TE0(m),TEs,R2s(m))).^2);
                
                % Parameter covariance matrices
                errOLS(:,:,m)=(D'*D)\(D'*M*D)/(D'*D); % unit weights
                errWLS(:,:,m)=pinv(D'/M*D); % optimal weights (M^-1)
            end
            
            % Check R2sEst and R2s are equal to within Nstd standard 
            % deviations of the theoretical error.
            % Nstd = 5 seems a good choice as it is fairly robust to 
            % outliers.
            Nstd=5;
            
            % R2*
            assertEqual(testCase,R2sEstOLS,R2s,'AbsTol',Nstd*reshape(sqrt(errOLS(1,1,:)),size(R2s)))
            assertEqual(testCase,R2sEstWLS,R2s,'AbsTol',Nstd*reshape(sqrt(errWLS(1,1,:)),size(R2s)))
            
            % Extrapolation to TE=0
            assertEqual(testCase,w_TE0_EstOLS{1},w_TE0,'AbsTol',w_TE0.*Nstd.*reshape(sqrt(errOLS(2,2,:)),size(w_TE0)))
            assertEqual(testCase,w_TE0_EstWLS{1},w_TE0,'AbsTol',w_TE0.*Nstd.*reshape(sqrt(errWLS(2,2,:)),size(w_TE0)))

        end
        
        function testMultipleContrastOLSvsWLS(testCase)
            
            dims=[100,10,10];
            R2s=100*rand(dims)+50; % in [50,150] / s
            
            % If max TE is too long, then estimation becomes inaccurate because of the
            % exponential decrease in SNR with TE
            TEs1=(2:2.5:20)*1e-3; % s
            
            w1_TE0=2000*rand(dims)+500;
            w1=hmri_test_utils.createDecaySignal(w1_TE0,TEs1,R2s);
                        
            % Add noise
            sigma1=1; % sigma1 = 1 => SNR = w1_TE0
            w1N=abs(w1+sigma1*randn(size(w1)));
            
            TEs2=TEs1(1:6); % s
            w2_TE0=2000*rand(dims)+500;
            w2=hmri_test_utils.createDecaySignal(w2_TE0,TEs2,R2s);
            
            % Add noise
            sigma2=1; % sigma2 = 1 => SNR = w2_TE0
            w2N=abs(w2+sigma2*randn(size(w2)));
            
            [R2sEstOLS,w_TE0_EstOLS]=hmri_calc_R2s([struct('data',w1N,'TE',TEs1),struct('data',w2N,'TE',TEs2)],'OLS');
            [R2sEstWLS,w_TE0_EstWLS]=hmri_calc_R2s([struct('data',w1N,'TE',TEs1),struct('data',w2N,'TE',TEs2)],'WLS1');
            
            % Compute theoretical covariance matrix of parameter estimates
            % https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation
            D=[TEs1(:),ones(length(TEs1),1),zeros(length(TEs1),1);
                TEs2(:),zeros(length(TEs2),1),ones(length(TEs2),1)];
            errOLS=zeros([3,3,dims]);
            errWLS=zeros([3,3,dims]);
            for m=1:prod(dims)
                % Covariance matrix of additive noise after log transform
                M=diag([sigma1./hmri_test_utils.createDecaySignal(w1_TE0(m),TEs1,R2s(m)),...
                    sigma2./hmri_test_utils.createDecaySignal(w2_TE0(m),TEs2,R2s(m))].^2);
                
                % Parameter covariance matrices
                errOLS(:,:,m)=(D'*D)\(D'*M*D)/(D'*D); % unit weights
                errWLS(:,:,m)=pinv(D'/M*D); % optimal weights (M^-1)
            end
            
            % Check R2sEst and R2s are equal to within Nstd standard 
            % deviations of the theoretical error.
            % Nstd = 5 seems a good choice as it is fairly robust to 
            % outliers.
            Nstd=5;
            
            % R2*
            assertEqual(testCase,R2sEstOLS,R2s,'AbsTol',Nstd*reshape(sqrt(errOLS(1,1,:)),size(R2s)))
            assertEqual(testCase,R2sEstWLS,R2s,'AbsTol',Nstd*reshape(sqrt(errWLS(1,1,:)),size(R2s)))
            
            % Extrapolation to TE=0 (w1)
            assertEqual(testCase,w_TE0_EstOLS{1},w1_TE0,'AbsTol',w1_TE0.*Nstd.*reshape(sqrt(errOLS(2,2,:)),size(w1_TE0)))
            assertEqual(testCase,w_TE0_EstWLS{1},w1_TE0,'AbsTol',w1_TE0.*Nstd.*reshape(sqrt(errWLS(2,2,:)),size(w1_TE0)))
            
            % Extrapolation to TE=0 (w2)
            assertEqual(testCase,w_TE0_EstOLS{2},w2_TE0,'AbsTol',Nstd*w2_TE0.*reshape(sqrt(errOLS(3,3,:)),size(w2_TE0)))
            assertEqual(testCase,w_TE0_EstWLS{2},w2_TE0,'AbsTol',Nstd*w2_TE0.*reshape(sqrt(errWLS(3,3,:)),size(w2_TE0)))

        end
    end
    
    %% Test Setup and Teardown Functions
    
    methods(TestMethodSetup)
        % These methods are run before each test
        function seedRandomNumberGenerator(testCase) %#ok<MANU>
            hmri_test_utils.seedRandomNumberGenerator;
        end
    end

    methods(TestClassSetup)
        % These methods are run when instantiating the class
        function enableRngLegacyWarnings(testCase) %#ok<MANU>
            % Make sure the use of legacy code that breaks random number
            % handling is going to generate a warning
            warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
            warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState') 
        end
    end
 
end
