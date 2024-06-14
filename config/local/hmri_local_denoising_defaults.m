% PURPOSE
% To set user-defined (site- or protocol-specific) defaults parameters
% which are used by the denoising module in the hMRI toolbox. Customized 
% processing parameters can be defined, overwriting defaults from 
% hmri_local_denoising defaults.
%__________________________________________________________________________
% Written by B. Ugurcan, 2024.

% init the global variable which carries the params
global hmri_def;

% Enter denoising default values as: hmri_def.denoising.(denoising-protocol).(denoising-parameter)

% The default values for lcpca denoising protocol
% all optional parameters turned off
hmri_def.denoising.lcpca_denoise.min_dimension =  0; % required-initialize here
hmri_def.denoising.lcpca_denoise.max_dimension = -1; % required-initialize here
hmri_def.denoising.lcpca_denoise.unwrap = false; % optional
hmri_def.denoising.lcpca_denoise.rescale_phs = false; % optional
hmri_def.denoising.lcpca_denoise.process_2d = false; % optional
hmri_def.denoising.lcpca_denoise.use_rmt = false; % optional
