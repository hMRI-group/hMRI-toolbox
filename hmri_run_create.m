function out = hmri_run_create(job)
%==========================================================================
% PURPOSE
% Calculation of multiparameter maps using B1 maps for B1 bias correction.
% If no B1 maps available, one can choose not to correct for B1 bias or
% apply UNICORT.
%==========================================================================

out.R1 = {};
out.R2s = {};
out.A = {};
out.MT = {};
out.T1w = {};
out.MTw = {};
out.PDw = {};

% loop over subjects in the main function, calling the local function for
% each subject:
for in=1:numel(job.subj)
    local_job.subj = job.subj(in);
    out_temp       = hmri_create_local(local_job);
    out.subj(in)   = out_temp.subj(1);
    out.R1{end+1}  = out.subj(in).R1{1};
    out.R2s{end+1} = out.subj(in).R2s{1};
    out.MT{end+1}  = out.subj(in).MT{1};
    out.A{end+1}   = out.subj(in).A{1};
    out.T1w{end+1} = out.subj(in).T1w{1};
    out.MTw{end+1} = out.subj(in).MTw{1};
    out.PDw{end+1} = out.subj(in).PDw{1};
end
end

%% =======================================================================%
% LOCAL SUBFUNCTION (PROCESSING FOR ONE SUBJET)
%=========================================================================%
function out_loc = hmri_create_local(job)

% determine output directory path
try 
    outpath = job.subj.output.outdir{1}; % case outdir
    if ~exist(outpath,'dir'); mkdir(outpath); end
catch  %#ok<CTCH>
    Pin = char(job.subj.raw_mpm.PD);
    outpath = fileparts(Pin(1,:)); % case indir
end
% save outpath as default for this job
hmri_get_defaults('outdir',outpath);

% Directory structure for results:
% <output directory>/Results
% <output directory>/Results/Supplementary
% <output directory>/B1mapCalc
% <output directory>/RFsensCalc
% <output directory>/MPMCalc
% The *Calc directories are deleted at the end of the Map Creation
% processing if hmri.cleanup defaults is set to true.
% If repeated runs, <output directory> is replaced by <output
% directory>/Run_xx to avoid overwriting previous outputs.

% define a directory for final results
% RESULTS contains the 4 final maps which are the essentials for the users
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
% SUPPLEMENTARY (within the Results directory) contains useful
% supplementary files (processing parameters and a few additional maps)
supplpath = fullfile(outpath, 'Results', 'Supplementary');
if ~exist(supplpath,'dir'); mkdir(supplpath); end

% define other (temporary) paths for processing data
b1path = fullfile(outpath, 'B1mapCalc');
if ~exist(b1path,'dir'); mkdir(b1path); end
rfsenspath = fullfile(outpath, 'RFsensCalc');
if ~exist(rfsenspath,'dir'); mkdir(rfsenspath); end
mpmpath = fullfile(outpath, 'MPMCalc');
if ~exist(mpmpath,'dir'); mkdir(mpmpath); end

% save all these paths in the job.subj structure
job.subj.path.b1path = b1path;
job.subj.path.b1respath = supplpath; % copy B1 maps to supplementary folder
job.subj.path.rfsenspath = rfsenspath;
job.subj.path.mpmpath = mpmpath;
job.subj.path.respath = respath;
job.subj.path.supplpath = supplpath;

% save log file location
job.subj.log.logfile = fullfile(supplpath, 'hMRI_map_creation_logfile.log');
job.subj.log.flags = struct('LogFile',struct('Enabled',true,'FileName','hMRI_map_creation_logfile.log','LogDir',supplpath), ...
    'PopUp',job.subj.popup,'ComWin',true);
flags = job.subj.log.flags;
flags.PopUp = false;
hmri_log(sprintf('\t============ CREATE hMRI MAPS MODULE - %s.m (%s) ============', mfilename, datetime('now')),flags);

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
spm_jsonwrite(fullfile(supplpath,'hMRI_map_creation_job_create_maps.json'),job,struct('indent','\t'));

% run B1 map calculation for B1 bias correction
P_trans = hmri_create_b1map(job.subj);

% check, if RF sensitivity profile was acquired and do the recalculation
% accordingly
if isfield(job.subj.sensitivity,'RF_once') || isfield(job.subj.sensitivity,'RF_per_contrast')
  job.subj = hmri_create_RFsens(job.subj);
end

% run hmri_create_MTProt to evaluate the parameter maps
job.subj.b1_trans_input = P_trans;
[fR1, fR2s, fMT, fA, PPDw, PT1w, PMTw]  = hmri_create_MTProt(job.subj);

% collect outputs:
out_loc.subj.R1  = {fR1};
out_loc.subj.R2s = {fR2s};
out_loc.subj.MT  = {fMT};
out_loc.subj.A   = {fA};
out_loc.subj.T1w = {PT1w};
out_loc.subj.MTw = {PMTw};
out_loc.subj.PDw = {PPDw};

% clean after if required
if hmri_get_defaults('cleanup')
    rmdir(job.subj.path.b1path,'s');
    rmdir(job.subj.path.rfsenspath,'s');
    rmdir(job.subj.path.mpmpath,'s');
end

f = fopen(fullfile(respath, '_finished_'), 'wb');
fclose(f);

hmri_log(sprintf('\t============ CREATE hMRI MAPS MODULE: completed (%s) ============', datetime('now')),flags);

end