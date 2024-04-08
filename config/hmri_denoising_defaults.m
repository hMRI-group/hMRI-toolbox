function hmri_denoising_defaults

%init the global variable which carries the params
global hmri_def;

%Enter denoising default values as: hmri_def.denoising.(denoising-protocol).(default-value)

%The default values for lcpca denoising protocol
%all optional parameters turned off
hmri_def.denoising.lcpca_denoise.mag_input = {}; %%required-null initialize here input with GUI
hmri_def.denoising.lcpca_denoise.phase_input = {}; %%optional-null initialize here input with GUI
hmri_def.denoising.lcpca_denoise.output_path=''; %%required-null initialize here input with GUI
hmri_def.denoising.lcpca_denoise.min_dimension= 0; %%required-initialize here 
hmri_def.denoising.lcpca_denoise.max_dimension = -1; %%required-initialize here 
hmri_def.denoising.lcpca_denoise.unwrap = false; %%optional
hmri_def.denoising.lcpca_denoise.rescale_phs = false; %%optional
hmri_def.denoising.lcpca_denoise.process_2d=false; %%optional
hmri_def.denoising.lcpca_denoise.use_rmt=false; %%optional
end