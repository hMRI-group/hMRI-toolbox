%% Main function to generate tests
function tests = hmri_calc_R2s_test
tests = functiontests(localfunctions);
end

%% Test Functions
function testSingleContrast3D(testCase)

% Test on 3D simulated data that the calculated R2* is within a defined
% tolerance given a single contrast input.

dims=[32,77,68];
R2s=100*rand(dims)+50; % in [50,150] / s

TEs=(2:2.5:20)*1e-3; % s
signal_TE0=2000*rand(dims)+500;

% Create signal decay given R2s and TEs:
signal=decaySignal(signal_TE0,TEs,R2s);

% Check that it works without outputting the extrapolated values
R2sEst=hmri_calc_R2s(struct('data',signal,'TE',TEs));

% Check every element is less than pre-defined threshold, here 1e-9 %TODO
% -rationale for threshold.
assertLessThan(testCase,abs(R2s-R2sEst),1e-9)

end


function testMultipleContrast3D(testCase)

% Test on 3D simulated data that the calculated R2* and contrast-specific
% intercepts are within a defined tolerance given multiple contrast input.

dims=[46,57,62];
R2s=100*rand(dims)+50; % in [50,150] / s

TEs1=(2:2.5:20)*1e-3; % s
signal1_TE0=2000*rand(dims)+500; % [500 2500]
signal1=decaySignal(signal1_TE0,TEs1,R2s);

% First four TEs for second contrast
TEs2=TEs1(1:4); % s
signal2_TE0=1000*rand(dims)+100; % [100 1100]
signal2=decaySignal(signal2_TE0,TEs2,R2s);

[R2sEst,extrapolated]=hmri_calc_R2s([struct('data',signal1,'TE',TEs1),struct('data',signal2,'TE',TEs2)]);

assertLessThan(testCase,abs(R2s-R2sEst),1e-9)
assertLessThan(testCase,abs(extrapolated{1}-signal1_TE0),1e-9)
assertLessThan(testCase,abs(extrapolated{2}-signal2_TE0),1e-9)

end

function testSingleContrast1D(testCase)

% Test on 1D simulated data that the calculated R2* and intercept
% are within a defined tolerance given single contrast input.

dims=[32,1];
R2s=100*rand(dims)+50; % in [50,150] / s

TEs=(2:2.5:20)*1e-3; % s
signal_TE0=2000*rand(dims)+500;
signal=decaySignal(signal_TE0,TEs,R2s);

[R2sEst,extrapolated]=hmri_calc_R2s(struct('data',signal,'TE',TEs));

assertLessThan(testCase,abs(R2s-R2sEst),1e-9)
assertLessThan(testCase,abs(extrapolated{1}-signal_TE0),1e-9)

end

function testZero2DInputs(testCase)

% Test on 2D simulated data that the calculated R2* is NaN for a single
% contrast input where the signal is actually zero.

dims=[32,56,8]; % last dim is echoes
signal=zeros(dims);
TEs=(2:2.5:20)*1e-3; % s

% Check that it fails when zeros are input
R2sEst=hmri_calc_R2s(struct('data',signal,'TE',TEs));

% Check every element is less than pre-defined threshold, here 1e-9 %TODO
% -rationale for threshold.
assertTrue(testCase,all(isnan(R2sEst(:))))

end

% function testNonStructInput(testCase)
% 
% import matlab.unittest.constraints.Throws;
% 
% TEs=(2:2.5:20)*1e-3; % s
% 
% % TO DO Check that an error is thrown if the input isn't a struct
% % testCase.assertThat(hmri_calc_R2s(zeros([32,56,8])), Throws('hmri:structError'));
%     
% end

function test1DNoise(testCase)

dims=[32,1];
R2s=100*rand(dims)+50; % in [50,150] / s

% If max TE is too long, then estimation becomes inaccurate because of the
% exponential decrease in SNR with TE
TEs=(2:2.5:20)*1e-3; % s
w_TE0=2000*rand(dims)+500;
w=decaySignal(w_TE0,TEs,R2s);

% SNR = w_TE0
wN=w+randn(size(w));

[R2sEst,extrapolated]=hmri_calc_R2s(struct('data',wN,'TE',TEs));

% Check relative error less than 5%
assertLessThan(testCase,abs(R2s-R2sEst)./R2s,0.05)
assertLessThan(testCase,abs(extrapolated{1}-w_TE0)./w_TE0,0.05)

end



function testMultipleContrastNoise(testCase)

dims=[32,77,68];
R2s=100*rand(dims)+50; % in [50,150] / s

% If max TE is too long, then estimation becomes inaccurate because of the
% exponential decrease in SNR with TE
TEs1=(2:2.5:20)*1e-3; % s

w1_TE0=2000*rand(dims)+500;
w1=decaySignal(w1_TE0,TEs1,R2s);
w1N=w1+randn(size(w1));

TEs2=TEs1(1:6); % s
w2_TE0=2000*rand(dims)+500;
w2=decaySignal(w2_TE0,TEs2,R2s);
w2N=w2+randn(size(w2));

[R2sEst,extrapolated]=hmri_calc_R2s([struct('data',w1N,'TE',TEs1),struct('data',w2N,'TE',TEs2)]);

% Check relative error less than 5%
assertLessThan(testCase,abs(R2s-R2sEst)./R2s,0.05)
assertLessThan(testCase,abs(extrapolated{1}-w1_TE0)./w1_TE0,0.05)
assertLessThan(testCase,abs(extrapolated{2}-w2_TE0)./w2_TE0,0.05)

end

function w_TEs=decaySignal(w_TE0,TEs,R2s)

dims=size(w_TE0);

% Account for 1D case
if (length(dims)==2)&&(dims(2)==1), dims=dims(1); end

TEs=reshape(TEs,[ones(1,length(dims)),length(TEs)]);

w_TEs=w_TE0.*exp(-R2s.*TEs);

end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end
