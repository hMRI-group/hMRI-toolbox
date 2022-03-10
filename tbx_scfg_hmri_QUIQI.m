function quiqi = tbx_scfg_hmri_QUIQI
% 
% 1./
% PURPOSE: Compute dictionaries of covariance matrices from the motion 
% degradation index.The dictionary  will be used by spm_reml.m or 
% spm_reml_sc.m to account for heteroscedasticity of the data.
% See " reference to the paper " for more details
% 
%
% METHODS: SPM.mat created when designing the model is modified to add in
% SPM.xVi.Vi the dictionary of covariance matrices. % 
% 
% 2./
% PURPOSE: Plot Residuals of the estimated model with respect to the Motion
% Degradation Index to evaluate the efficiency of the weighting
% 
%
% METHODS: Spatial variance of the residuals of each subject is plotted
% with respect to their corresponfing MDI value
%
%
%_______________________________________________________________________
% Nadege Corbin
% 2021.03.30
% Centre de Resonance Magnetique des Systemes Biologiques, Bordeaux, France

quiqi        = cfg_choice;
quiqi.tag     = 'quiqi';
quiqi.name    = 'QUIQI';
quiqi.help    = {'Compute dictionaries of covariance matrices from the motion',... 
' degradation index.The disctionary  will be used by spm_reml.m or',... 
' spm_reml_sc.m to account for heteroscedasticity of the data.'
    }'; 
quiqi.values  = {tbx_scfg_hmri_QUIQI_build tbx_scfg_hmri_QUIQI_check};


end
