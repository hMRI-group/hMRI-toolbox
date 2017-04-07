function out = hmri_run_mpr_b0_b1(job)
% Calculation of B1 mapping data 3D EPI spin echo (SE) and stimulated (STE) 
% echo images (see Jiru and Klose MRM 2006).
% Corresponding scanning protocol/sequence: al_B1mapping_v2a and
% al_B1mapping_v2b
% Input: 11 pairs of (SE, STE) images for B1 map calculation and 3 images
% for B0 map calculation.
% This macro calls the functions hmri_B1Map_unwarp and hmri_B1Map_process for
% correction of image distortions, padding and smoothing of the images.
% Output:
%     - distorted B1 (B1map_----) and error (SDmap_----) maps
%     - undistorted B1 (uB1map_----) and error (uSDmap_----) maps
%     - undistorted, masked and padded B1 maps (muB1map_---------)
%     - undistorted, masked, padded and smoothed B1 maps (smuB1map_---------) i.e. FULLY PROCESSED
% At each voxel, this macro selects the 5 pairs of (SE,STE image) (out of
% 11) with maximum signal amplitude in the SE images.
% The sum of square image of all SE images is created (SumOfSq) and
% undistorted (uSumOfSq) for coregistration of the B1 map to an anatomical dataset
% former hmri_B1map_v2.m

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
    out_temp       = hmri_mpr_b0_b1_local(local_job);
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

function out_loc = hmri_mpr_b0_b1_local(job)

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

% check, if RF sensitivity profile was acquired and do the recalculation
% accordingly
if ~isfield(job.subj.sensitivity,'RF_none')
  job.subj = hmri_RFsens(job.subj);
end

% initialize file names for map creation
P_mtw    = char(job.subj.raw_mpm.MT);
P_pdw    = char(job.subj.raw_mpm.PD);
P_t1w    = char(job.subj.raw_mpm.T1);

P_receiv = [];

% run hmri_MTProt to evaluate the parameter maps
[fR1, fR2s, fMT, fA, PPDw, PT1w]  = hmri_MTProt(P_mtw, P_pdw, P_t1w, P_trans, P_receiv);

out_loc.subj.R1  = {fullfile(outpath,spm_str_manip(fR1,'t'))};
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