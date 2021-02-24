%% Main function to generate tests
function tests = hmri_calc_R2s_test
tests = functiontests(localfunctions);
end

%% Test Functions
function testSingleContrast(testCase)

dims=[32,77,68];
R2s=100*rand(dims)+50; % in [50,150] / s

TEs=(2:2.5:20)*1e-3; % s
w=2000*rand(dims)+500;
w=decaySignal(w,TEs,R2s);

% Check that it works without outputting the extrapolated values
R2sEst=hmri_calc_R2s(struct('data',w,'TE',TEs));

assertLessThan(testCase,abs(R2s-R2sEst),1e-9)

end

function test1D(testCase)

dims=[32,1];
R2s=100*rand(dims)+50; % in [50,150] / s

TEs=(2:2.5:20)*1e-3; % s
w_TE0=2000*rand(dims)+500;
w=decaySignal(w_TE0,TEs,R2s);

[R2sEst,extrapolated]=hmri_calc_R2s(struct('data',w,'TE',TEs));

assertLessThan(testCase,abs(R2s-R2sEst),1e-9)
assertLessThan(testCase,abs(extrapolated{1}-w_TE0),1e-9)

end

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

function testMultipleContrast(testCase)

dims=[32,77,68];
R2s=100*rand(dims)+50; % in [50,150] / s

TEs1=(2:2.5:20)*1e-3; % s

w1_TE0=2000*rand(dims)+500;
w1=decaySignal(w1_TE0,TEs1,R2s);

TEs2=TEs1(1:4); % s
w2_TE0=2000*rand(dims)+500;
w2=decaySignal(w2_TE0,TEs2,R2s);

[R2sEst,extrapolated]=hmri_calc_R2s([struct('data',w1,'TE',TEs1),struct('data',w2,'TE',TEs2)]);

assertLessThan(testCase,abs(R2s-R2sEst),1e-9)
assertLessThan(testCase,abs(extrapolated{1}-w1_TE0),1e-9)
assertLessThan(testCase,abs(extrapolated{2}-w2_TE0),1e-9)

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
