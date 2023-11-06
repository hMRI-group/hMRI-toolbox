% PURPOSE
%   These are the suggested defaults for correcting MTsat for B1-inhomogeneity in the 
%   postmortem brain protocol from Lipp, et al. (MRM 2023, https://doi.org/10.1002/mrm.29524)
% 
% CHANGES 
%   Set correction model to the linear "lipp" model

% Global hmri_def variable used across the whole toolbox
global hmri_def

%--------------------------------------------------------------------------
% Which model to use for B1-correction of MTsat
%--------------------------------------------------------------------------
hmri_def.MTsatB1CorrectionModel = 'lipp'; % 'helms' or 'lipp'
hmri_def.MTsatB1CorrectionLippC = 1.2; % value for 700° MT pulse in postmortem brain protocol in Lipp, et al. (MRM 2023)
