function hmri_denoising_defaults
% Sets the defaults parameters which are used by the denoising module of the
% hMRI toolbox.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% DON'T MODIFY THIS FILE, IT CONTAINS THE REFERENCE DEFAULTS PARAMETERS.
% Please refer to hMRI-Toolbox\config\local\hmri_local_denoising_defaults.m
% to customise defaults parameters.  
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% FORMAT hmri_denoising_defaults
%__________________________________________________________________________
%
% THIS FILE SHOULD NOT BE MODIFIED. 
% To customize the hMRI-Toolbox defaults parameters so they match your own
% site- or protocol-specific setup, please refer to the defaults files in
% hMRI-Toolbox\config\local. In particular, use 
% "hmri_local_denoising_defaults.m".Make a copy with meaningful name, 
% modify as desired and select it in the hMRI toolbox denoising GUI.
%
% The structure and content of this file are largely inspired by the
% equivalent file in SPM.
%
% DOCUMENTATION
% A brief description of each parameter is provided together with
% guidelines and recommendations to modify these parameters. With few
% exceptions, parameters should ONLY be MODIFIED and customized BY ADVANCED
% USERS, having a good knowledge of the underlying algorithms and
% implementation. Please refer to the documentation in the github WIKI for
% more details.  
%__________________________________________________________________________
% Written by B. Ugurcan, 2024.

% init the global variable which carries the params
global hmri_def;

% Enter denoising default values as: hmri_def.denoising.(denoising-protocol).(denoising-parameter)

% The default values for lcpca denoising protocol
hmri_def.denoising.lcpca_denoise.min_dimension =  0; % minimum number of kept PCA components
hmri_def.denoising.lcpca_denoise.max_dimension = -1; % maximum number of kept PCA components (-1 for all components)
hmri_def.denoising.lcpca_denoise.unwrap = true; % whether to unwrap the phase data or keep it as is
hmri_def.denoising.lcpca_denoise.rescale_phs = true; % whether to rescale the phase data or keep it as is, assuming radians
hmri_def.denoising.lcpca_denoise.process_2d = false; % whether to denoise in 2D, for instance when acquiring a thin slab of data
hmri_def.denoising.lcpca_denoise.use_rmt = false; % whether to use random matrix theory rather than noise fitting to estimate the noise threshold

end