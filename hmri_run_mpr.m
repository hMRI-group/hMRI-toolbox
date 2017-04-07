function out = hmri_run_mpr(job)
% Calculation of multiparameter maps using B1 maps for B1 bias correction.
% If no B1 maps available, one can choose not to correct for B1 bias or
% apply UNICORT.

job = hmri_process_data_spec(job);

out.R1 = {};
out.R2s = {};
out.A = {};
out.MT = {};
out.T1w = {};

% loop over subjects in the main function, calling the local function for
% each subject:
for in=1:numel(job.subj)
    local_job.subj = job.subj(in);
    out_temp       = hmri_mpr_local(local_job);
    out.subj(in)   = out_temp.subj(1);
    out.R1{end+1}  = out.subj(in).R1{1};
    out.R2s{end+1} = out.subj(in).R2s{1};
    out.MT{end+1}  = out.subj(in).MT{1};
    out.A{end+1}   = out.subj(in).A{1};
    out.T1w{end+1} = out.subj(in).T1w{1};
end
end

% ========================================================================
%% SUBFUNCTION
% ========================================================================

function out_loc = hmri_mpr_local(job)

% determine output directory path
try 
    outpath = job.subj.output.outdir{1}; % case outdir
catch 
    Pin = char(job.subj.raw_mpm.MT);
    outpath = fileparts(Pin(1,:)); % case outdir
end
% save outpath as default for this job
hmri_get_defaults('outdir',outpath);

% run B1 map calculation for B1 bias correction
P_trans = hmri_run_b1map(job.subj);

% initialize file names for map creation
P_mtw    = char(job.subj.raw_mpm.MT);
P_pdw    = char(job.subj.raw_mpm.PD);
P_t1w    = char(job.subj.raw_mpm.T1);

P_receiv = [];

% run hmri_MTProt to evaluate the parameter maps
[fR1, fR2s, fMT, fA, PPDw, PT1w]  = hmri_MTProt(P_mtw, P_pdw, P_t1w, P_trans, P_receiv);

% apply UNICORT if required, and collect output:
if strcmp(job.subj.b1_type,'UNICORT')
    out_unicort = hmri_run_unicort(PPDw, fR1);
    out_loc.subj.R1  = {fullfile(outpath,spm_str_manip(out_unicort.R1u,'t'))};
else
    out_loc.subj.R1  = {fullfile(outpath,spm_str_manip(fR1,'t'))};
end
out_loc.subj.R2s = {fullfile(outpath,spm_str_manip(fR2s,'t'))};
out_loc.subj.MT  = {fullfile(outpath,spm_str_manip(fMT,'t'))};
out_loc.subj.A   = {fullfile(outpath,spm_str_manip(fA,'t'))};
out_loc.subj.T1w = {fullfile(outpath,spm_str_manip(PT1w,'t'))};

% save processing params (hmri defaults) and job for the current subject:
hmri_def = hmri_get_defaults;
save(fullfile(outpath, [spm_file(P_mtw(1,:),'basename') '_create_maps_hmridef.mat']),'hmri_def');
save(fullfile(outpath, [spm_file(P_mtw(1,:),'basename') '_create_maps_job.mat']),'job');

f = fopen(fullfile(outpath, '_finished_'), 'wb');
fclose(f);

end