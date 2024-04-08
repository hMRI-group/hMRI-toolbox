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
hmri_log(sprintf('\t============ DENOISING MODULE - %s.m (%s) ============', mfilename, datestr(now)),flags);

if newrespath
    hmri_log(sprintf(['WARNING: existing results from previous run(s) were found, \n' ...
        'the output directory has been modified. It is now:\n%s\n'],outpath),job.subj.log.flags);
else
    hmri_log(sprintf('INFO: the output directory is:\n%s\n',outpath),flags);
end

%check basic requirements and error if basic reqirements fail
switch dntype{1}
    case 'lcpca_denoise'
        check_params.mag_input = cellstr(char(spm_file(job.subj.denoisingtype.(dntype{1}).mag_input,'number','')));
        check_params.phase_input = cellstr(char(spm_file(job.subj.denoisingtype.(dntype{1}).phase_input,'number','')));
        check_params.phase_bool = any(~cellfun(@isempty, check_params.phase_input));
        check_params.mag_bool = any(~cellfun(@isempty, check_params.mag_input));

        %Issue an error and abort in cases of non-existent magnitude image data and
        %non-equal number of non-empty phase and magnitude images

        if ~check_params.mag_bool || ~isfield(job.subj.denoisingtype.(dntype{1}),'mag_input')
        hmri_log('No magnitude images were entered, aborting! Please check your input data and try again!', flags);
        error('Error: No magnitude images were entered, aborting! Please check your input data and try again!')
        end

        if check_params.phase_bool && (length(check_params.mag_input) ~= length(check_params.phase_input))
        hmri_log('The number of phase and magnitude images are different, please check your input data and try again!', flags);
        error("Error: The number of phase and magnitude images are different, please check your input data and try again!")
        end
end
 
% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
job.SPMver = sprintf('%s (%s)', v, r);

% save original job
spm_jsonwrite(fullfile(supplpath,'hMRI_denoising_job.json'),job,struct('indent','\t'));


end