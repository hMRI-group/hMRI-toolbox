function out = hmri_run_create_B1(job)
%==========================================================================
% PURPOSE
% Calculation of B1 maps for B1 bias correction.
%==========================================================================

% loop over subjects in the main function, calling the local function for
% each subject:
for in=1:numel(job.subj)
    local_job.subj = job.subj(in);
    out_temp       = hmri_create_local(local_job);
    out.subj(in)   = out_temp.subj(1);
end
end

%% =======================================================================%
% LOCAL SUBFUNCTION (PROCESSING FOR ONE SUBJECT)
%=========================================================================%
function out_loc = hmri_create_local(job)

% determine output directory path
b1type=fields(job.subj.b1_type);
assert(~isempty(job.subj.b1_type.(b1type{1}).b1input),'Data to make B1 map must be provided!')
try 
    outpath = job.subj.output.outdir{1}; % case outdir
    if ~exist(outpath,'dir'); mkdir(outpath); end
catch  %#ok<CTCH>
    Pin = char(job.subj.b1_type.(b1type{1}).b1input);
    outpath = fileparts(Pin(1,:)); % case indir
end
% save outpath as default for this job
hmri_get_defaults('outdir',outpath);

% Directory structure for results:
% <output directory>/Results
% <output directory>/Results/Supplementary
% <output directory>/B1mapCalc
% If repeated runs, <output directory> is replaced by <output
% directory>/Run_xx to avoid overwriting previous outputs.

% define a directory for final results
% RESULTS contains the final maps which are the essentials for the users
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

% define other paths for processing data
b1path = fullfile(outpath, 'B1mapCalc');
if ~exist(b1path,'dir'); mkdir(b1path); end

% save all these paths in the job.subj structure
job.subj.path.b1path = b1path;
job.subj.path.b1respath = respath; % copy B1 maps to results folder
job.subj.path.respath = respath;
job.subj.path.supplpath = supplpath;

% save log file location
job.subj.log.logfile = fullfile(supplpath, 'hMRI_B1_map_creation_logfile.log');
job.subj.log.flags = struct('LogFile',struct('Enabled',true,'FileName','hMRI_B1_map_creation_logfile.log','LogDir',supplpath), ...
    'PopUp',job.subj.popup,'ComWin',true);
flags = job.subj.log.flags;
flags.PopUp = false;
hmri_log(sprintf('\t============ CREATE B1 MAP MODULE - %s.m (%s) ============', mfilename, datetime('now')),flags);

if newrespath
    hmri_log(sprintf(['WARNING: existing results from previous run(s) were found, \n' ...
        'the output directory has been modified. It is now:\n%s\n'],outpath),job.subj.log.flags);
else
    hmri_log(sprintf('INFO: the output directory is:\n%s\n',outpath),flags);
end

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
job.SPMver = sprintf('%s (%s)', v, r);

% save original job (before it gets modified by RFsens)
spm_jsonwrite(fullfile(supplpath,'hMRI_map_creation_job_create_B1_map.json'),job,struct('indent','\t'));

% run B1 map calculation for B1 bias correction
P_trans = cellstr(hmri_create_b1map(job.subj));

% collect outputs:
out_loc.subj.B1ref=P_trans(1);
out_loc.subj.B1map=P_trans(2);

out_loc.subj.B1  = P_trans;

% clean after if required
if hmri_get_defaults('cleanup')
    rmdir(job.subj.path.b1path,'s');
end

f = fopen(fullfile(respath, '_finished_'), 'wb');
fclose(f);

hmri_log(sprintf('\t============ CREATE B1 MAP MODULE: completed (%s) ============', datetime('now')),flags);

end
