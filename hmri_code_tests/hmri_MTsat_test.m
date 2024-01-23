%% Main function to generate tests
function tests = hmri_MTsat_test
tests = functiontests(localfunctions);
end

%% Test Functions

function typical3TprotocolTest(testCase)

MTw.TR=23.7e-3; % s
MTw.fa=deg2rad(6); % rad

MTw=addDataFields(MTw,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map,testCase.TestData.MTsat/100);

MTsatEst=hmri_calc_MTsat(MTw, testCase.TestData.PD, testCase.TestData.R1);

assertEqual(testCase,MTsatEst,testCase.TestData.MTsat,'AbsTol',testCase.TestData.tol,'Estimated MTsat has large error!')

end

function helmsB1correctionTest(testCase)

MTsat=100*testCase.TestData.MTsat;
B1=testCase.TestData.B1map;
C=0.4;

MTsatB1bias=MTsat.*(1-C.*B1)/(1-C);

MTsatEst=hmri_correct_MTsat(MTsatB1bias, B1, 'helms', C);

assertEqual(testCase,MTsatEst,MTsat,'AbsTol',testCase.TestData.tol,'Estimated MTsat has large error!')

end

function lippB1correctionTest(testCase)

MTsat=testCase.TestData.MTsat;
B1=testCase.TestData.B1map;
C=1.2;

MTsatB1bias=MTsat.*(1+C.*(B1-1));

MTsatEst=hmri_correct_MTsat(MTsatB1bias, B1, 'lipp', C);

assertEqual(testCase,MTsatEst,MTsat,'AbsTol',testCase.TestData.tol,'Estimated MTsat has large error!')

end

%% subfunctions
function w=addDataFields(w,PD,R1,B1,delta)

w.data=bsxfun(@times,PD,hmri_test_utils.HelmsMTwApprox(w.fa*B1,w.TR,R1,delta));
w.B1=B1;

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

% Reset random seed for reproducibility
hmri_test_utils.seedRandomNumberGenerator;

% Define reasonable parameter set
dims=[20,50,30];

R1min=0.5; % / s
R1max=2; % / s
testCase.TestData.R1=R1min+(R1max-R1min)*rand([dims,1]); % / s

PDmin=500;
PDmax=2000;
testCase.TestData.PD=PDmin+(PDmax-PDmin)*rand([dims,1]);

B1min=0.4;
B1max=1.6;
testCase.TestData.B1map=B1min+(B1max-B1min)*rand([dims,1]);

MTsatMin=0.1;
MTsatMax=2;
testCase.TestData.MTsat=MTsatMin+(MTsatMax-MTsatMin)*rand([dims,1]);

testCase.TestData.tol=1e-3*(MTsatMin+MTsatMax/2);

end

function teardown(testCase)  % do not change function name
clear MTw MTsatEst MTsat MTsatB1bias B1 cfg_check_assignin
end