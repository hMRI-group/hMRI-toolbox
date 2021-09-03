%% Main function to generate tests
function tests = hmri_calc_A_test
tests = functiontests(localfunctions);
end

%% Test Functions
function typical7TprotocolTest(testCase)
PDw.TR=28.5e-3; % s
PDw.fa=deg2rad(5); % rad
T1w.TR=28.5e-3; % s
T1w.fa=deg2rad(26); % rad

small_angle_approx=false;

PDw=addDataFields(PDw,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map);
T1w=addDataFields(T1w,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map);

PDest=hmri_calc_A(PDw,T1w,small_angle_approx);

verifyLessThan(testCase,abs(PDest(:)-testCase.TestData.PD(:)),testCase.TestData.tol,'Estimated PD has large error!')

end

function typical3TprotocolTest(testCase)
PDw.TR=25e-3; % s
PDw.fa=deg2rad(6); % rad
T1w.TR=25e-3; % s
T1w.fa=deg2rad(21); % rad

small_angle_approx=false;

PDw=addDataFields(PDw,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map);
T1w=addDataFields(T1w,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map);

PDest=hmri_calc_A(PDw,T1w,small_angle_approx);

verifyLessThan(testCase,abs(PDest(:)-testCase.TestData.PD(:)),testCase.TestData.tol,'Estimated PD has large error!')

end

function compareToOldMethodTest(testCase)

[newerr,olderr]=compareToOldMethod(testCase,false);
verifyLessThan(testCase,abs(newerr)-abs(olderr),testCase.TestData.tol,'Small angle approximation gives significantly smaller residuals than the estimation without this approximation in some cases!')

[newerr,olderr]=compareToOldMethod(testCase,true);
verifyLessThan(testCase,abs(newerr)-abs(olderr),testCase.TestData.tol,'Old calculation method gives significantly smaller residuals than the new one in some cases!')

end

%% subfunctions
% Ernst equation
function S=ernst(alpha,TR,R1)

S=sin(alpha).*(1-exp(-TR.*R1))./(1-cos(alpha).*exp(-TR.*R1));

end

function w=addDataFields(w,PD,R1,B1)

w.data=reshape(bsxfun(@times,PD,ernst(w.fa*B1(:),w.TR,R1)),size(B1));
w.B1=B1;

end

function [newerr,olderr]=compareToOldMethod(testCase,small_angle_approx)

% Protocol from Weiskopf, et al. (2013).
PDw.TR=23.7e-3; % s
PDw.fa=deg2rad(6); % rad
T1w.TR=18.7e-3; % s
T1w.fa=deg2rad(20); % rad

PDw=addDataFields(PDw,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map);
T1w=addDataFields(T1w,testCase.TestData.PD,testCase.TestData.R1,testCase.TestData.B1map);

PDest=hmri_calc_A(PDw,T1w,small_angle_approx);

% Old implementation of R1 calculation in hMRI toolbox
PDsa=zeros(size(PDest));
for p = 1:size(PDw.data,3)
    
    PDwSlice = PDw.data(:,:,p);
    T1wSlice = T1w.data(:,:,p);
    
    f_T = testCase.TestData.B1map(:,:,p); 
    
    % Transmit bias corrected quantitative T1 values
    % correct T1 for transmit bias f_T with fa_true = f_T * fa_nom
    % T1corr = T1 / f_T / f_T
    T1 = ((((PDwSlice / PDw.fa) - (T1wSlice / T1w.fa)+eps) ./ ...
        max((T1wSlice * T1w.fa / 2 / T1w.TR) - (PDwSlice * PDw.fa / 2 / PDw.TR),eps))./f_T.^2);
    
    PDsa(:,:,p) = T1 .* (T1wSlice .*(T1w.fa*f_T) / 2 / T1w.TR) + (T1wSlice ./ (T1w.fa*f_T));
    
end

newerr=PDest(:)-testCase.TestData.PD(:);
olderr=PDsa(:)-testCase.TestData.PD(:);

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

% Define reasonable parameter set
dims=[20,50,30];

R1min=0.5; % / s
R1max=2; % / s
testCase.TestData.R1=R1min+(R1max-R1min)*rand(prod(dims),1); % / s

PDmin=500;
PDmax=2000;
testCase.TestData.PD=PDmin+(PDmax-PDmin)*rand(prod(dims),1);

B1min=0.4;
B1max=1.6;
testCase.TestData.B1map=B1min+(B1max-B1min)*rand([dims,1]);

testCase.TestData.tol=1e-3*(PDmin+PDmax)/2;

end

function teardown(testCase)  % do not change function name
clear PDw T1w PDest errinerr PDsa
end