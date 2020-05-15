function hmri_corr_imperf_spoil(job)
%==========================================================================
% PURPOSE
% Compute correction factors for Imperfect spoiling
% As described in (Preibisch & Deichmann, MRM,2009)
% Numerical simulations are performed using the EPG formalism 
% EPG dunctions come from the toolbox https://github.com/mriphysics/EPG-X,
% from Shaihan Malik, KCL
% The parameters used for the simulation and the resulting correction
% factors are written in a json file in the output folder selected by the 
% user.
%==========================================================================


%% ***********************************************%%
% 1./ Numerical simulations with EPG
%*************************************************%%
%%
% Get sequence parameters
FA=job.seq_params.FA_deg; % (FA1 FA2) [deg]
TR=job.seq_params.TR_ms;  % [ms]
Phi0=job.seq_params.Phi0_deg; % [deg]
B1range=job.seq_params.B1range; % normalized to 1
Gdur=job.seq_params.Gdur_ms; % [ms]
Gamp=job.seq_params.Gamp_mt_per_m; % [[mt/m]

if length(Gdur)~= length (Gamp)
    error('The vector of Gradient Durations and Gradient Amplitude must have the size')
end

%% Get tissue parameters
T1range=job.tissue_params.T1range_ms; %[ms]
T2=job.tissue_params.T2_ms; % [ms]
D=job.tissue_params.D_um2_per_ms; % [um^2/ms]

%% MRI constants
gamma = 42.58e6*2*pi; % rad.Hz/T

%% Build studture "diff" to account for diffusion effect
diff = struct;
diff.D = D*1e-9;
diff.G=Gamp;
diff.tau=Gdur;

%% Run EPG simulation
S = zeros([2 length(T1range) length(B1range)]);
for T1val = 1 : length(T1range) % loop over T1 values
    
    T1 = T1range(T1val);
    npulse = floor(15*T1/TR);   % To ensure steady state signal
    nstart = ceil(10*T1/TR);    % To specify a range over which to average (there is no average NC)
    
    for B1val = 1 : length(B1range)  % loop over B1+ values
        B1eff = B1range(B1val);
        
        %%% make train of flips and phases
        phi_train = RF_phase_cycle(npulse,Phi0); % phase of the RF pulses
        alpha_train1 = d2r(FA(1).*B1eff)*ones([1 npulse]); % flip angles of the PDw acquisitions
        alpha_train2 = d2r(FA(2).*B1eff)*ones([1 npulse]); % flip angles of the T1w acquisitions
        
        % Calculate signals via EPG:
        
        %PDw
        F0 = EPG_GRE(alpha_train1,phi_train,TR, T1, T2, 'diff', diff);
        S(1,T1val,B1val) = (abs(F0(end)));
        %T1w
        F0 = EPG_GRE(alpha_train2,phi_train,TR, T1, T2, 'diff', diff);
        S(2,T1val,B1val) = (abs(F0(end)));
        
    end
end



%% ***********************************************%%
% 2./ Fitting T1=A(B1eff)+B(B1eff)*T1app
%*************************************************%%

ABcoeff = zeros(2, length(B1range));
T1app = zeros(length(B1range), length(T1range));
for B1val = 1 : length(B1range)
    
    B1eff = B1range(B1val);
    
    % Calculate T1app, accounting for B1+
    T1app(B1val,:) = squeeze(2.*(S(1,:,B1val)./d2r(FA(1).*B1eff)-S(2,:,B1val)./d2r(FA(2).*B1eff))./(S(2,:,B1val).*d2r(FA(2).*B1eff)./TR-S(1,:,B1val).*d2r(FA(1).*B1eff)./TR));
    
    % build matrix X with column of ones and column of T1app
    X = ones([length(T1range) 2]);
    X(:,2) = T1app(B1val,:);
    ABcoeff(:, B1val) = pinv(X)*T1range';
    
end

%% *********************************************************%%
% 3./ Fitting A=P(B1eff) and B=P(B1eff) with 2nd degree polynom
%***********************************************************%%

polyCoeffA = polyfit(B1range, ABcoeff(1,:), 2);
polyCoeffB = polyfit(B1range, ABcoeff(2,:), 2);

%% *********************************************************%%
% 4./ Compute RMSE on T1app and T1
%***********************************************************%%
T1corr = repmat(polyval(polyCoeffA, B1range)',[1 length(T1range)]) + repmat(polyval(polyCoeffB, B1range)',[1 length(T1range)]).*T1app;
T1_Corr_Err = (T1corr - repmat(T1range, [length(B1range) 1]))./repmat(T1range, [length(B1range) 1])*100;
T1_App_Err = (T1app - repmat(T1range, [length(B1range) 1]))./repmat(T1range, [length(B1range) 1])*100;

RMSE_Corr = sqrt(mean(T1_Corr_Err(:).^2));
RMSE_App = sqrt(mean(T1_App_Err(:).^2));

%% *********************************************************%%
% 5./ Write parameters and correction factors in a json file in the selected
% output directory
%***********************************************************%%
Results.Input=job;
Results.Output.P2_a = round(polyCoeffA,3);
Results.Output.P2_b = round(polyCoeffB,3);
Results.Output.RMSE_percent.T1app=round(RMSE_App,3);
Results.Output.RMSE_percent.T1corr=round(RMSE_Corr,3);

Results.ToCopy{1}=['hmri_def.MPMacq_set.names{1} = ' job.prot_name ];
Results.ToCopy{end+1}=['hmri_def.MPMacq_set.tags{1}  = ' strrep(job.prot_name,' ','')];
Results.ToCopy{end+1}=['hmri_def.MPMacq_set.vals{1}  = [' num2str([TR TR FA]) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.tag = ' strrep(job.prot_name,' ','') ];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_a = [' num2str(round(polyCoeffA,3)) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_b = [' num2str(round(polyCoeffB,3)) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.enabled = hmri_def.imperfectSpoilCorr.enabled;'];

results_filename=fullfile(job.outdir,strrep(job.prot_name,' ',''));

spm_jsonwrite(results_filename{1},Results,struct('indent','\t'));
end