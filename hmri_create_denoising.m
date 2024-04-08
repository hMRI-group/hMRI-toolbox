function varargout = hmri_create_denoising(job)

% retrieve effective acquisition & processing parameters
jobsubj = job.subj;
denoisedout = get_denoising_params(jobsubj);
protocolfield = fieldnames(jobsubj.denoisingtype);
denoising_protocol = protocolfield{1};
end