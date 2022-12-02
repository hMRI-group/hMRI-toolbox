%% Main function to generate tests
function tests = hmri_make_errormaps_test
tests = functiontests(localfunctions);
end

function hmri_make_dR1_test(testCase)
htu = hmri_test_utils;

% Define reasonable parameter set
dims=[20,50,30];

% Make sure that always the same random number is produced
htu.seedRandomNumberGenerator;

R1min=0.5; % / s
R1max=2; % / s
R1=R1min+(R1max-R1min)*rand((dims)); % / s

PDmin=100;
PDmax=1000;
PD=PDmin+(PDmax-PDmin)*rand((dims));

B1min=0.4;
B1max=1.6;
B1map=B1min+(B1max-B1min)*rand([dims,1]);

tol=1e-2;

%% Typical 7T protocol
PDw.TR=28.5e-3; % s
PDw.fa=deg2rad(5); % rad
PDw.B1=B1map;
T1w.TR=28.5e-3; % s
T1w.fa=deg2rad(26); % rad
T1w.B1=B1map;

PDw.data=reshape(bsxfun(@times,PD,htu.ernst(PDw.fa*B1map,PDw.TR,R1)),dims);
T1w.data=reshape(bsxfun(@times,PD,htu.ernst(T1w.fa*B1map,T1w.TR,R1)),dims);

sigma=0.1;

nrep = 10;

errorinPD=sigma*randn([dims nrep]);
PDwErr=PDw;
PDwErr.data=PDw.data+errorinPD;

errorinT1=sigma*randn([dims nrep]);
T1wErr=T1w;
T1wErr.data=T1w.data+errorinT1;

R1est=hmri_calc_R1(PDwErr,T1wErr,B1map);

errorinR1=hmri_make_dR1(PDwErr.data,T1wErr.data,errorinPD,errorinT1,PDwErr.fa,T1wErr.fa,PDwErr.TR,T1wErr.TR,B1map,R1est);

merrorinR1 = mean(errorinR1,4);

rms_error = std(R1-R1est,[],4);

% figure
% hold on
% histogram(rms_error)
% histogram(merrorinR1)

assert(rms((rms_error(:)-merrorinR1(:))./R1(:))<tol)

end


function hmri_make_dPD_test(testCase)
htu = hmri_test_utils;

% Define reasonable parameter set
dims=[20,50,30];

% make sure that always the same random number is produced
htu.seedRandomNumberGenerator;

R1min=0.5; % / s
R1max=2; % / s
R1=R1min+(R1max-R1min)*rand((dims)); % / s

PDmin=100;
PDmax=1000;
PD=PDmin+(PDmax-PDmin)*rand((dims));

B1min=0.4;
B1max=1.6;
B1map=B1min+(B1max-B1min)*rand([dims,1]);

tol=1e-2;

%% Typical 7T protocol
PDw.TR=28.5e-3; % s
PDw.fa=deg2rad(5); % rad
PDw.B1=B1map;
T1w.TR=28.5e-3; % s
T1w.fa=deg2rad(26); % rad
T1w.B1=B1map;

PDw.data=reshape(bsxfun(@times,PD,htu.ernst(PDw.fa*B1map,PDw.TR,R1)),dims);
T1w.data=reshape(bsxfun(@times,PD,htu.ernst(T1w.fa*B1map,T1w.TR,R1)),dims);

sigma=0.1;

nrep = 10;

errorinPD=sigma*randn([dims nrep]);
PDwErr=PDw;
PDwErr.data=PDw.data+errorinPD;

errorinT1=sigma*randn([dims nrep]);
T1wErr=T1w;
T1wErr.data=T1w.data+errorinT1;

Aest=hmri_calc_A(PDwErr,T1wErr,B1map);

modelerrorinPD=hmri_make_dPD(PDwErr.data,T1wErr.data,errorinPD,errorinT1,PDwErr.fa,T1wErr.fa,PDwErr.TR,T1wErr.TR,B1map,Aest);

mmodelerrorinPD = mean(modelerrorinPD,4);

rms_error = std(PD-Aest,[],4);

% figure
% hold on
% histogram(rms_error)
% histogram(mmodelerrorinPD)
% rms((rms_error(:)-mmodelerrorinPD(:))./PD(:))

assert(rms((rms_error(:)-mmodelerrorinPD(:))./PD(:))<tol)

end