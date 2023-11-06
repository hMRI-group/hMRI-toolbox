% PURPOSE
%   These are the suggested defaults for correcting MTsat for B1-inhomogeneity when using the 7T 
%   protocol from Olsson, et al. (MRM 2021, https://doi.org/10.1002/mrm.28899)
% 
% CHANGES 
%   Set correction model to the linear "helms" model and set parameter to the optimal value from 
%   Olsson, et al. (MRM 2021)

% Global hmri_def variable used across the whole toolbox
global hmri_def

%--------------------------------------------------------------------------
% Which model to use for B1-correction of MTsat
%--------------------------------------------------------------------------
hmri_def.MTsatB1CorrectionModel = 'helms'; % 'helms' or 'lipp'
hmri_def.MTsatB1CorrectionHelmsC = 0.34; % value for 7T protocol from Olsson, et al. (MRM 2021)
