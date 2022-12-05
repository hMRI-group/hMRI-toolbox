function hmri_local_defaults_errormaps
% PURPOSE
% To set user-defined (site- or protocol-specific) defaults parameters
% which are used by the hMRI toolbox. Customized processing parameters can
% be defined, overwriting defaults from hmri_defaults. Acquisition
% protocols can be specified here as a fallback solution when no metadata
% are available. Note that the use of metadata is strongly recommended. 
%
% RECOMMENDATIONS
% Parameters defined in this file are identical, initially, to the ones
% defined in hmri_defaults.m. It is recommended, when modifying this file,
% to remove all unchanged entries and save the file with a meaningful name.
% This will help you identifying the appropriate defaults to be used for
% each protocol, and will improve the readability of the file by pointing
% to the modified parameters only.
%
% WARNING
% Modification of the defaults parameters may impair the integrity of the
% toolbox, leading to unexpected behaviour. ONLY RECOMMENDED FOR ADVANCED
% USERS - i.e. who have a good knowledge of the underlying algorithms and
% implementation. The SAME SET OF DEFAULT PARAMETERS must be used to
% process uniformly all the data from a given study. 
%
% HOW DOES IT WORK?
% The modified defaults file can be selected using the "Configure toolbox"
% branch of the hMRI-Toolbox. For customization of B1 processing
% parameters, type "help hmri_b1_standard_defaults.m". 
%
% DOCUMENTATION
% A brief description of each parameter is provided together with
% guidelines and recommendations to modify these parameters. With few
% exceptions, parameters should ONLY be MODIFIED and customized BY ADVANCED
% USERS, having a good knowledge of the underlying algorithms and
% implementation. 
% Please refer to the documentation in the github WIKI for more details. 
%__________________________________________________________________________
% This is an example defaults file for the error map calculation.
% S.Mohammadi 12/07/2022

% Global hmri_def variable used across the whole toolbox
global hmri_def

%==========================================================================
% Generation of error maps
%==========================================================================
% errormaps
hmri_def.errormaps  = true;

% lower bounds on R1, PD and MTsat to avoid dividing by small numbers
% in parameter-wise SNR estimates of R1/dR1, PD/dPD and MTsat/dMTsat
hmri_def.qMRI_maps_thresh.dR1 = 1e-4;
hmri_def.qMRI_maps_thresh.dPD = 1e-2;
hmri_def.qMRI_maps_thresh.dMT = 1e-4;

% upper bounds on SNR maps dR1/R1, dPD/PD and dMTsat/MTsat
hmri_def.qMRI_maps_thresh.SNR_R1= 1e3;
hmri_def.qMRI_maps_thresh.SNR_PD= 1e3;
hmri_def.qMRI_maps_thresh.SNR_MT= 1e3;

% weighted combination - see Mohammadi et al., NeuroImage, 2022 for details
% might have to be adjusted for different dataset
hmri_def.wcombparams.kt         = 10;   % This parameter is relevant and might have to be adjusted per protocol (defined in percent). 
% Examples can be found in Mohammadi et al., NeuroImage, 2022, Supplementary Materials: S1: Efficiency of robust combination and the Fermi function
hmri_def.wcombparams.am         = true; % True: arithmetic mean is written out, in addition to robust combination
hmri_def.wcombparams.errormaps  = true; % True: robust combination of error maps is generated and written out
% The following parameters are for experts only.
hmri_def.wcombparams.res        = -4;   % Resampling factor, determines interpolation method (details can be found in spm_slice_vol.m)
hmri_def.wcombparams.smthk      = 0;    % If larger than zero, the error maps will be spatially smoothed with the a Gaussian kernel of FWHM = smthk
hmri_def.wcombparams.dim        = 3;    % Defines dimention along which slices along which weightes will be determined

end
