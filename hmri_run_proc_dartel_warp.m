function out = hmri_run_proc_dartel_warp(job)
% Wrapper function for the DARTEL processing part of the hMRI toolbox.
% This relies on the standard SPM functions but handles some extra
% settings, e.g. the selection of the output folder.
% For this, things run as usual with SPM functions then data are moved
% around where they should be.

% Fix the job structure
job_spm = job;

% Define output folder for SPM function
output = struct('outDir', [], 'option', 'same');
if isfield(job.output,'outdir') % -> everything in the same
    output.outDir = job.output.outdir{1};
    output.option = 'allin';
elseif isfield(job.output,'outdir_ps') % -> per suject organization
    output.outDir = job.output.outdir_ps{1};
    output.option = 'subjspec';
end
job_spm.output = output;

% Call SPM's 'run' function
out = spm_dartel_warp(job_spm);

end
