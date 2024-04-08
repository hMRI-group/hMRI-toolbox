function out = hmri_run_denoise(job)

%Get denoising protocol
dntype=fields(job.subj.denoisingtype);

%Case indir versus outdir
try 
    outpath = job.subj.output.outdir{1}; % case outdir
    if ~exist(outpath,'dir'); mkdir(outpath); end
catch  
    Pin = char(job.subj.denoisingtype.(dntype{1}).mag_input);
    outpath = fileparts(Pin(1,:)); % case indir
end


end