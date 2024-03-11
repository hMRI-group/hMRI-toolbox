function hmri_corr_imperf_spoil(job)
%==========================================================================
% PURPOSE
% Compute coefficients to correct for the effect of imperfect spoiling on
% T1 estimation as described in (Preibisch & Deichmann, MRM 2009).
%
% Numerical simulations are performed using the EPG formalism described in
% Malik et al., MRM 2017 and available here:
% https://github.com/mriphysics/EPG-X
%
% The parameters used for the simulation and the resulting correction
% factors are written to a JSON file in the output folder specified by the
% user.
%==========================================================================

hmri_log(sprintf('\t--- Calculating Imperfect Spoiling Correction Coefficients ---'));

if isfield(job.params,'params_diff')
    hmri_log(sprintf('\t......with diffusion spoiling'));
    diff_flag = true;

    seq_params    = job.params.params_diff.seq_params_diff;
    tissue_params = job.params.params_diff.tissue_params_diff;

    %% Build structure "diff" to account for diffusion spoiling
    diff     = struct;
    diff.D   = tissue_params.D_um2_per_ms*1e-9;
    diff.G   = seq_params.Gamp_mT_per_m;
    diff.tau = seq_params.Gdur_ms;

    assert(length(diff.G) == length(diff.tau), 'The vectors of gradient durations and amplitudes must have the same length!')
    

elseif isfield(job.params,'params_nodiff')
    hmri_log(sprintf('\t......without diffusion spoiling'));
    diff_flag = false;

    seq_params    = job.params.params_nodiff.seq_params;
    tissue_params = job.params.params_nodiff.tissue_params;

else
    error('something went wrong!')
end

%% ***********************************************%%
% 1./ Numerical simulations with EPG
%*************************************************%%
%%
% Get sequence parameters
FA      = seq_params.FA_deg;              % Flip angles [deg]
TR      = seq_params.TR_ms;               % [ms]
Phi0    = seq_params.Phi0_deg;            % [deg]
B1range = seq_params.B1range_percent/100; % convert such that 100% = 1

%% Get tissue parameters
T1range = tissue_params.T1range_ms;     %[ms]
T2range = tissue_params.T2range_ms;     % [ms]


%% Run EPG simulation
nT1 = length(T1range);
nT2 = length(T2range);
nB1 = length(B1range);
S1  = zeros([nT1 nT2 nB1]);
S2  = zeros([nT1 nT2 nB1]);
hmri_log(sprintf('\t-------- Simulating signals'));
for T1val = 1 : nT1 % loop over T1 values, can use parfor for speed

    T1 = T1range(T1val);
    npulse = floor(15*T1/min(TR));   % To ensure steady state signal

    for T2val = 1 : nT2
        T2 = T2range(T2val);

        for B1val = 1 : nB1  % loop over B1+ values
            B1eff = B1range(B1val);

            %%% make train of flip angles and their phases
            phi_train = RF_phase_cycle(npulse,Phi0); % phase of the RF pulses
            alpha_train1 = d2r(FA(1).*B1eff)*ones([1 npulse]); % flip angles of the PDw acquisitions
            alpha_train2 = d2r(FA(2).*B1eff)*ones([1 npulse]); % flip angles of the T1w acquisitions

            % Calculate signals via EPG:
            if diff_flag
                %PDw
                F0 = EPG_GRE(alpha_train1, phi_train, TR(1), T1, T2, 'diff', diff);
                S1(T1val,T2val,B1val) = (abs(F0(end)));
                %T1w
                F0 = EPG_GRE(alpha_train2, phi_train, TR(2), T1, T2, 'diff', diff);
                S2(T1val,T2val,B1val) = (abs(F0(end)));
            else
                %PDw
                F0 = EPG_GRE(alpha_train1, phi_train, TR(1), T1, T2);
                S1(T1val,T2val,B1val) = (abs(F0(end)));
                %T1w
                F0 = EPG_GRE(alpha_train2, phi_train, TR(2), T1, T2);
                S2(T1val,T2val,B1val) = (abs(F0(end)));
            end

        end
    end
end



%% ***********************************************%%
% 2./ Fitting T1=A(B1eff)+B(B1eff)*T1app
%*************************************************%%
hmri_log(sprintf('\t-------- Determining Coefficients'));
ABcoeff = zeros(2, nB1);
T1app = zeros(nB1, nT1, nT2);
for B1val = 1 : nB1

    B1eff = B1range(B1val);

    % Calculate T1app, accounting for B1+
    T1app(B1val,:,:) = 1./hmri_calc_R1(...
        struct('data',S1(:,:,B1val),'fa',d2r(FA(1)),'TR',TR(1),'B1',B1eff),...
        struct('data',S2(:,:,B1val),'fa',d2r(FA(2)),'TR',TR(2),'B1',B1eff),...
        job.small_angle_approx);

    % build matrix X with column of ones and column of T1app
    X = ones([nT1*nT2 2]);
    X(:,2) = T1app(B1val,:);
    ABcoeff(:, B1val) = pinv(X)*repmat(T1range, [1 nT2]).';

end

%% *********************************************************%%
% 3./ Fitting A=P(B1eff) and B=P(B1eff) with 2nd degree polynom
%***********************************************************%%
polyCoeffA = polyfit(B1range, ABcoeff(1,:), 2);
polyCoeffB = polyfit(B1range, ABcoeff(2,:), 2);


%% *********************************************************%%
% 4./ Compute RMSE on T1app and T1
%***********************************************************%%
T1app = T1app(:,:);
T1corr = repmat(polyval(polyCoeffA, B1range).',[1 nT1*nT2])+ repmat(polyval(polyCoeffB, B1range).',[1 nT1*nT2]).*T1app;
T1_Corr_Err = (T1corr - repmat(T1range, [nB1 nT2]))./repmat(T1range, [nB1 nT2])*100;
T1_App_Err = (T1app - repmat(T1range, [nB1 nT2]))./repmat(T1range, [nB1 nT2])*100;


RMSE_Corr = sqrt(mean(T1_Corr_Err(:).^2));
RMSE_App = sqrt(mean(T1_App_Err(:).^2));


%% *********************************************************%%
% 5./ Write parameters and correction factors in a json file
% in the selected output directory
%***********************************************************%%
hmri_log(sprintf('\t-------- Writing results\n'));
Results.Input = job;
Results.Output.P2_a = round(polyCoeffA,4);
Results.Output.P2_b = round(polyCoeffB,4);
Results.Output.small_angle_approx = job.small_angle_approx;
Results.Output.RMSE_percent.T1app=round(RMSE_App,3);
Results.Output.RMSE_percent.T1corr=round(RMSE_Corr,3);

Results.ToCopy{1}    =['hmri_def.MPMacq_set.names{NN} = ''' job.prot_name ''';' ];
Results.ToCopy{end+1}=['hmri_def.MPMacq_set.tags{NN}  = ''' strrep(job.prot_name,' ','') ''';'];
Results.ToCopy{end+1}=['hmri_def.MPMacq_set.vals{NN}  = [' num2str([TR FA]) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.tag = ''' strrep(job.prot_name,' ','') ''';' ];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_a = [' num2str(round(polyCoeffA,4)) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_b = [' num2str(round(polyCoeffB,4)) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.small_angle_approx = ' mat2str(job.small_angle_approx) ';'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.enabled = hmri_def.imperfectSpoilCorr.enabled;'];

results_filename = fullfile(job.outdir,[strrep(job.prot_name,' ',''),'.json']);

spm_jsonwrite(results_filename{1},Results,struct('indent','\t'));

end