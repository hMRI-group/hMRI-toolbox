function hmri_MTsat_Lipp_model_500deg_defaults
% PURPOSE
%   These are the suggested defaults for correcting MTsat for B1-inhomogeneity when using the standard Siemens MT-module
% 
% CHANGES 
%   Set correction model to the linear "Lipp" model and scaled the coefficient from a referance MT flip angle of 700° to 500° 

% Global hmri_def variable used across the whole toolbox
global hmri_def

%--------------------------------------------------------------------------
% Which model to use for B1-correction of MTsat
%--------------------------------------------------------------------------
hmri_def.MTsatB1CorrectionModel = 'lipp'; % 'helms' or 'lipp'
hmri_def.MTsatB1CorrectionLippC = 1.3; % value for 500° scaled from the 700° value in Lipp, et al. (MRM 2023) using Eq. (S.3) in the supplementary material to that paper

end
