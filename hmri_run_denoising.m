function out = hmri_run_denoising(job)

% loop over subjects in the main function, calling the local function for
% each subject:
for in=1:numel(job.subj)
    out.subj(in) = hmri_denoising_local(job.subj(in));
end

end

function out = hmri_denoising_local(jobsubj)

% concatenate all contrasts into jobsubj.mag_img and phase_img
contrasts = {'pdw','t1w','mtw'};
jobsubj.mag_img = {};
jobsubj.phase_img = {};
for c = 1:length(contrasts)
    con = contrasts{c};
    jobsubj.mag_img = [jobsubj.mag_img; jobsubj.(con).mag_img];
    jobsubj.phase_img = [jobsubj.phase_img; jobsubj.(con).phase_img];
end

% Case indir versus outdir
try
    outpath = jobsubj.output.outdir{1}; % case outdir
    if ~exist(outpath,'dir'); mkdir(outpath); end
catch
    Pin = char(jobsubj.mag_img);
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

% save all these paths in the jobsubj structure
jobsubj.path.dnrespath = respath;
jobsubj.path.respath = respath;
jobsubj.path.supplpath = supplpath;

% save log file location
jobsubj.log.logfile = fullfile(supplpath, 'hMRI_denoising_logfile.log');
jobsubj.log.flags = struct('LogFile',struct('Enabled',true,'FileName','hMRI_denoising_logfile.log','LogDir',supplpath), ...
    'PopUp',jobsubj.popup,'ComWin',true);
flags = jobsubj.log.flags;
flags.PopUp = false;
hmri_log(sprintf('\t============ DENOISING MODULE - %s.m (%s) ============', mfilename, datetime('now')),flags);

if newrespath
    hmri_log(sprintf(['WARNING: existing results from previous run(s) were found, \n' ...
        'the output directory has been modified. It is now:\n%s\n'],outpath),jobsubj.log.flags);
else
    hmri_log(sprintf('INFO: the output directory is:\n%s\n',outpath),flags);
end

% check basic requirements and error if basic requirements fail
check_params.mag_img = cellstr(char(spm_file(jobsubj.mag_img,'number','')));
check_params.mag_bool = any(~cellfun(@isempty, check_params.mag_img));
check_params.phase_img = cellstr(char(spm_file(jobsubj.phase_img,'number','')));
check_params.phase_bool = any(~cellfun(@isempty, check_params.phase_img));

% Issue an error and abort in cases of non-existent magnitude image data and
% non-equal number of non-empty phase and magnitude images
if ~check_params.mag_bool || ~isfield(jobsubj,'mag_img')
    msg = 'No magnitude images were entered, aborting! Please check your input data and try again!';
    hmri_log(msg, flags);
    error(msg)
end

if check_params.phase_bool && (length(check_params.mag_img) ~= length(check_params.phase_img))
    msg = 'The number of phase and magnitude images are different, please check your input data and try again!';
    hmri_log(msg, flags);
    error(msg)
end

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
jobsubj.SPMver = sprintf('%s (%s)', v, r);

% save original job to supplementary directory
spm_jsonwrite(fullfile(supplpath,'hMRI_denoising_job.json'),jobsubj,struct('indent','\t'));

% run the denoising and collect the results
hmri_log(sprintf('\t============ %s - %s.m (%s) ============',"APPLYING DENOISING", mfilename, datetime('now')),flags);

% run the denoising
[output_mag, output_phase] = hmri_denoising(jobsubj);

% assign output dependencies to denoised output images of separate contrasts
iMag = 1;
iPhase = 1;
for c = 1:length(contrasts)
    con = contrasts{c};

    nMag = length(jobsubj.(con).mag_img);
    idxstr = ['DenoisedMagnitude' con];
    out.(idxstr) = output_mag(iMag:iMag+nMag-1);
    iMag = iMag + nMag;
    
    nPhase = length(jobsubj.(con).phase_img);
    idxstr = ['DenoisedPhase' con];
    out.(idxstr) = output_phase(iPhase:iPhase+nPhase-1);
    iPhase = iPhase + nPhase;
end

hmri_log(sprintf('\t============ DENOISING MODULE: completed (%s) ============', datetime('now')),flags);
end