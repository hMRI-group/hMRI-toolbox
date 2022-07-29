function hmri_DIB_defaults
% PURPOSE
% These are the recommended defaults for the example dataset 
% (Callaghan, et al. 2019. Example dataset for the hMRI toolbox. Data in 
% Brief 25, 104132 https://doi.org/10.1016/j.dib.2019.104132).
% 
% CHANGES 
% 1) define centre and scanner
% 2) make the segmentation less aggressive (because the defacing part of the 
%    data anonymisation procedure causes edges near the brain which interact
%    with segmentation) 
% 3) turn on imperfect spoiling correction.

% Global hmri_def variable used across the whole toolbox
global hmri_def

% Specify the research centre & scanner. Not mandatory.
hmri_def.centre = 'fil' ;
hmri_def.scanner = 'prisma' ;

%==========================================================================
% Default parameters for segmentation
% ADVANCED USERS ONLY!
% hmri_def.segment is effectively the job to be handed to spm_preproc_run
% By default, parameters are set to
% - create tissue class images (c*) in the native space of the source image
%   (tissue(*).native = [1 0]) for tissue classes 1-5
% - save both BiasField and BiasCorrected volume (channel.write = [1 1])
% - recommended values from SPM12 (October 2017)
%==========================================================================
% Make clean up less aggressive due to defacing causing edges near the edge
% of the brain
hmri_def.segment.warp.mrf = 0.5;
hmri_def.segment.warp.cleanup = 1;

%--------------------------------------------------------------------------
% MPM acquisition parameters and RF spoiling correction parameters
%--------------------------------------------------------------------------
% Given the possible confusion and resulting mistake (imperfect spoiling
% correction applied to the wrong sequence), when TR and FA values match
% one of the listed cases below, the option is disabled by default. 
% When enabling the imperfect spoiling correction, make sure the
% coefficients retrieved in the list below are definitely calculated for
% the protocol used!
hmri_def.imperfectSpoilCorr.enabled = true;

end
