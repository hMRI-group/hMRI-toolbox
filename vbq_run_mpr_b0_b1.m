function out = vbq_run_mpr_b0_b1(job)
% Calculation of B1 mapping data 3D EPI spin echo (SE) and stimulated (STE) 
% echo images (see Jiru and Klose MRM 2006).
% Corresponding scanning protocol/sequence: al_B1mapping_v2a and
% al_B1mapping_v2b
% Input: 11 pairs of (SE, STE) images for B1 map calculation and 3 images
% for B0 map calculation.
% This macro calls the functions vbq_B1Map_unwarp and vbq_B1Map_process for
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
% former vbq_B1map_v2.m

% $Id$

job = vbq_process_data_spec(job);

out.R1 = {};
out.R2s = {};
out.A = {};
out.MT = {};
out.T1w = {};

% loop over subjects in the main function, calling the local function for
% each subject:
for in=1:numel(job.subj)
    local_job.subj = job.subj(in);
    out_temp       = vbq_mpr_b0_b1_local(local_job);
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

function out_loc = vbq_mpr_b0_b1_local(job)

% run B1 map calculation for B1 bias correction
P_trans = vbq_run_b1map(job.subj);

% initialize file names for map creation
P_mtw    = char(job.subj.raw_mpm.MT);
P_pdw    = char(job.subj.raw_mpm.PD);
P_t1w    = char(job.subj.raw_mpm.T1);

P_receiv = [];

% determine output directory path
% CASE INDIR
outpath = fileparts(P_mtw(1,:)); 
% CASE OUTDIR
if isfield(job.subj.output,'outdir')
    if ~strcmp(outpath, job.subj.output.outdir{1})
         vbq_get_defaults('outdir',job.subj.output.outdir{1});
         outpath = vbq_get_defaults('outdir');
    end
end

% run vbq_MTProt to evaluate the parameter maps
[fR1, fR2s, fMT, fA, PPDw, PT1w]  = vbq_MTProt(P_mtw, P_pdw, P_t1w, P_trans, P_receiv);

% % determine output directory path
% % CASE INDIR
% outpath = fileparts(fR1); 
% % CASE OUTDIR
% if isfield(job.subj.output,'outdir')
%     if ~strcmp(outpath, job.subj.output.outdir{1})
%         % MFC: Only move files if a different directory is chosen - can't
%         % move a file to itself...
%         outpath = job.subj.output.outdir{1};
%         movefile(fR1,outpath);
%         movefile(fR2s,outpath);
%         movefile(fMT,outpath);
%         movefile(fA,outpath);
%         movefile(PPDw,outpath);
%         movefile(PT1w,outpath);
%     end
% end

out_loc.subj.R1  = {fullfile(outpath,spm_str_manip(fR1,'t'))};
out_loc.subj.R2s = {fullfile(outpath,spm_str_manip(fR2s,'t'))};
out_loc.subj.MT  = {fullfile(outpath,spm_str_manip(fMT,'t'))};
out_loc.subj.A   = {fullfile(outpath,spm_str_manip(fA,'t'))};
out_loc.subj.T1w = {fullfile(outpath,spm_str_manip(PT1w,'t'))};

% save processing params (vbq defaults) and job for the current subject:
vbq_def = vbq_get_defaults; %#ok<*NASGU>
save(fullfile(outpath, [spm_file(P_mtw(1,:),'basename') '_create_maps_vbqdef.mat']),'vbq_def');
save(fullfile(outpath, [spm_file(P_mtw(1,:),'basename') '_create_maps_job.mat']),'job');

f = fopen(fullfile(outpath, '_finished_'), 'wb');
fclose(f);

end