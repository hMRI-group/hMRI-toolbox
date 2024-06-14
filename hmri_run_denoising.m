function out = hmri_run_denoising(job)

% loop over subjects in the main function, calling the local function for
% each subject:
for in=1:numel(job.subj)
    out.subj(in) = hmri_denoising_local(job.subj(in));
end

end

function out = hmri_denoising_local(job)

% Get denoising protocol
dntype=fieldnames(job.subj.denoisingtype);

% Case indir versus outdir
try
    outpath = job.subj.output.outdir{1}; % case outdir
    if ~exist(outpath,'dir'); mkdir(outpath); end
catch
    Pin = char(job.subj.denoisingtype.(dntype{1}).mag_input);
    outpath = fileparts(Pin(1,:)); % case indir
end

% save outpath as default for this job
hmri_get_defaults('outdir',outpath);

% Directory structure for results:
% <output directory>/Results
% <output directory>/Results/Supplementary
% If repeated runs, <output directory> is replaced by <output
% directory>/Run_xx to avoid overwriting previous outputs.

% define a directory for final results
% RESULTS contains the final denoised maps
% Supplementary containes additonal info and log
respath = fullfile(outpath, 'Results');
newrespath = false;
if exist(respath,'dir')
    index = 1;
    tmpoutpath = outpath;
    while exist(tmpoutpath,'dir')
        index = index + 1;
        tmpoutpath = fullfile(outpath,sprintf('Run_%0.2d',index));
    end
    outpath = tmpoutpath;
    mkdir(outpath);
    respath = fullfile(outpath, 'Results');
    newrespath = true;
end
if ~exist(respath,'dir'); mkdir(respath); end
supplpath = fullfile(outpath, 'Results', 'Supplementary');
if ~exist(supplpath,'dir'); mkdir(supplpath); end

% save all these paths in the job.subj structure
job.subj.path.dnrespath = respath;
job.subj.path.respath = respath;
job.subj.path.supplpath = supplpath;

% save log file location
job.subj.log.logfile = fullfile(supplpath, 'hMRI_denoising_logfile.log');
job.subj.log.flags = struct('LogFile',struct('Enabled',true,'FileName','hMRI_denoising_logfile.log','LogDir',supplpath), ...
    'PopUp',job.subj.popup,'ComWin',true);
flags = job.subj.log.flags;
flags.PopUp = false;
hmri_log(sprintf('\t============ DENOISING MODULE - %s.m (%s) ============', mfilename, datetime('now')),flags);

if newrespath
    hmri_log(sprintf(['WARNING: existing results from previous run(s) were found, \n' ...
        'the output directory has been modified. It is now:\n%s\n'],outpath),job.subj.log.flags);
else
    hmri_log(sprintf('INFO: the output directory is:\n%s\n',outpath),flags);
end

% check basic requirements and error if basic requirements fail
check_params.mag_input = cellstr(char(spm_file(job.subj.denoisingtype.(dntype{1}).mag_input,'number','')));
check_params.phase_input = cellstr(char(spm_file(job.subj.denoisingtype.(dntype{1}).phase_input,'number','')));
check_params.phase_bool = any(~cellfun(@isempty, check_params.phase_input));
check_params.mag_bool = any(~cellfun(@isempty, check_params.mag_input));

% Issue an error and abort in cases of non-existent magnitude image data and
% non-equal number of non-empty phase and magnitude images
if ~check_params.mag_bool || ~isfield(job.subj.denoisingtype.(dntype{1}),'mag_input')
    msg = 'No magnitude images were entered, aborting! Please check your input data and try again!';
    hmri_log(msg, flags);
    error(msg)
end

if check_params.phase_bool && (length(check_params.mag_input) ~= length(check_params.phase_input))
    msg = 'The number of phase and magnitude images are different, please check your input data and try again!';
    hmri_log(msg, flags);
    error(msg)
end

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
job.SPMver = sprintf('%s (%s)', v, r);

% save original job to supplementary directory
spm_jsonwrite(fullfile(supplpath,'hMRI_denoising_job.json'),job,struct('indent','\t'));

% run the denoising and collect the results
hmri_log(sprintf('\t============ %s - %s.m (%s) ============',"APPLYING DENOISING", mfilename, datetime('now')),flags);

% concatenate all contrasts into jobsubj.denoisingtype.(denoising_protocol).mag_img and phase_img
contrasts = {'mtw','pdw','t1w'};
jobsubj.denoisingtype.(dntype{1}).mag_img = {};
jobsubj.denoisingtype.(dntype{1}).phase_img = {};
for c = 1:length(contrasts)
    con = contrasts{c};

    jobsubj.denoisingtype.(dntype{1}).mag_img = 
        [jobsubj.denoisingtype.(dntype{1}).mag_img;
        jobsubj.denoisingtype.(dntype{1}).(con).mag_img];

    jobsubj.denoisingtype.(dntype{1}).phase_img = 
        [jobsubj.denoisingtype.(dntype{1}).phase_img;
        jobsubj.denoisingtype.(dntype{1}).(con).phase_img];
end

% run the denoising
[output_mag, output_phase] = hmri_denoising(job);

% assign output dependencies to denoised output images of separate contrasts
iMag = 1;
iPhase = 1;
for c = 1:length(contrasts)
    con = contrasts{c};

    nMag = length(jobsubj.denoisingtype.(dntype{1}).(con).mag_img);
    idxstr = ['DenoisedMagnitude' con];
    out.subj.(idxstr) = output_mag(iMag:iMag+nMag);
    iMag = iMag + nMag;
    
    nPhase = length(jobsubj.denoisingtype.(dntype{1}).(con).phase_img);
    idxstr = ['DenoisedPhase' con];
    out.subj.(idxstr) = output_phase(iPhase:iPhase+nPhase);
    iPhase = iPhase + nPhase;
end

hmri_log(sprintf('\t============ DENOISING MODULE: completed (%s) ============', datetime('now')),flags);
end