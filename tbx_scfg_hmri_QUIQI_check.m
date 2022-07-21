function quiqi_check=tbx_scfg_hmri_QUIQI_check
% 
% PURPOSE: Plot Residuals of the estimated model with respect to the Motion
% Degradation Index to evaluate the efficiency of the weighting
% 
%
% METHODS: Spatial variance of the residuals of each subject is plotted
% with respect to their corresponfing MDI value
%
%_______________________________________________________________________
% Nadege Corbin
% 2021.03.30
% Centre de Resonance Magnetique des Systemes Biologiques, Bordeaux, France
% ======================================================================

% ---------------------------------------------------------------------
% SPM.mat file
% ---------------------------------------------------------------------
spm_mat_file         = cfg_files;
spm_mat_file.tag     = 'spm_mat_file';
spm_mat_file.name    = 'SPM.mat file';
spm_mat_file.help    = {'Select the SPM.mat file containing the design of the model.'};
spm_mat_file.ufilter = '^SPM.mat$';
spm_mat_file.num     = [1 1];

pow        = cfg_entry;
pow.tag     = 'power';
pow.name    = 'Fit power ';
pow.val     = {[4]};
pow.strtype = 'e';
pow.num     = [1 1];
pow.help    = {'Specify the power of  the polynomial fit.'};


quiqi_check        = cfg_exbranch;
quiqi_check.tag     = 'quiqi_check';
quiqi_check.name    = 'QUIQI CHECK';
quiqi_check.val     = { spm_mat_file pow};
quiqi_check.help    = {'Spatial variance of the residuals of each subject is plotted ',...
 'with respect to their corresponfing MDI value'};
quiqi_check.prog    = @hmri_quiqi_check;
