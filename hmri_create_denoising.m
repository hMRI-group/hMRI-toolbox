function varargout = hmri_create_denoising(job)

% retrieve effective acquisition & processing parameters
jobsubj = job.subj;
denoisedout = get_denoising_params(jobsubj);
protocolfield = fieldnames(jobsubj.denoisingtype);
denoising_protocol = protocolfield{1};

%execute the chosen denoising method and define output
switch denoising_protocol
    case 'lcpca_denoise'
[output_mag, output_phase] = hmri_calc_lcpcadenoise(denoisedout);
varargout{1} = output_mag;
varargout{2} = output_phase;
end

%=========================================================================%
% Write/Arrange denoising parameters tp the denoising_params.
%=========================================================================%
function denoising_params = get_denoising_params(jobsubj)

% get denoising method
dntypename = fieldnames(jobsubj.denoisingtype);
denoising_protocol = dntypename{1};

% init defaults filename and optional-defaults bool 
deffnam = '';
custom_def = false;
end


end