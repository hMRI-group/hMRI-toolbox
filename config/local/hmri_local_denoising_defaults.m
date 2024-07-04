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
hmri_def.denoising.lcpca_denoise.min_dimension =  0; % minimum number of kept PCA components
hmri_def.denoising.lcpca_denoise.max_dimension = -1; % maximum number of kept PCA components (-1 for all components)
hmri_def.denoising.lcpca_denoise.unwrap = true; % whether to unwrap the phase data or keep it as is
hmri_def.denoising.lcpca_denoise.rescale_phs = true; % whether to rescale the phase data or keep it as is, assuming radians
hmri_def.denoising.lcpca_denoise.process_2d = false; % whether to denoise in 2D, for instance when acquiring a thin slab of data
hmri_def.denoising.lcpca_denoise.use_rmt = false; % whether to use random matrix theory rather than noise fitting to estimate the noise threshold
