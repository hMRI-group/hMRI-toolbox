function out = hmri_run_proc_pipeline(job)
% Deal with the preprocessing pipelines. There are 2 options
% 1/ US+Smooth
% 2/ US+Dartel+Smooth
%
% Input include only some parametric maps, the structural maps (for
% segmentation), the required smoothing and which pipeline to use. All
% other options are hard-coded!
% By default, these pipelines focus only on the first 2 tissue classes,
% i.e. GM and WM only.
% 
% For more flexibility, individual modules can be combined. :-)
% 
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

% 1/ Setup the smoothing and execute hmri_run_proc_US
%----------------------------------------------------
% Get the proc_US job structure, with all the defaults
proc_us = tbx_scfg_hmri_proc_US;
[~, job_US] = harvest(proc_us, proc_us, false, false);

% Fill in the data now: parametric maps & structurals for segmentation
if job.pipe_c==1 % US+smooth -> warp the parametric maps here.
    job_US.many_sdatas.vols_pm = job.vols_pm;
end
job_US.many_sdatas.channel.vols = job.s_vols;
% job_US.many_sdatas.rstruct.s_vols = job.s_vols;

% Only spit out CSF in native space (no Dartel imported) & warped (no modulation)
job_US.tissue(3).native = [1 0];
job_US.tissue(3).warped = [1 0];

% Run the *_proc_US
fprintf('\nhMRI-pipeline: running the US module.\n')
out_US = hmri_run_proc_US(job_US);
% where the output structure 'out_US'
% .tiss : struct-array with subfields
%           .c and .rc, for the native and Dartel imported
%           .wc and .mwc, for the warped and modulated
%         tissue class images
% .maps : struct-array with subfields 'wvols_pm' for the warped parametric
%         maps
% .def  : cell-array with the deformations for each subject.

% 2/ Proceed with dartel (only of requested)
%-------------------------------------------
% including template create and warping into MNI space

if job.pipe_c == 2
    % DARTEL processing: 
    proc_Dartel = tbx_scfg_hmri_proc_Dartel;
    % a) warping with template creation
    proc_Dwarp = proc_Dartel.values{1};
    [~, job_Dwarp] = harvest(proc_Dwarp, proc_Dwarp, false, false);
    
    % Fill in the filenames for rc1 and rc2, i.e. *only* GM and WM
    for ii=1:2
        job_Dwarp.images{ii} = spm_file(out_US.tiss(ii).rc,'number',1);
    end
    % Run the Dartel-warp job
    fprintf('\nhMRI-pipeline: running the Darte-warp module.\n')
    out_Dwarp = spm_dartel_template(job_Dwarp);
    
    % b) normalize to  MNI
    proc_Dnorm = proc_Dartel.values{3};
    [~, job_Dnorm] = harvest(proc_Dnorm, proc_Dnorm, false, false);
    % get last tempalte
    job_Dnorm.template = out_Dwarp.template(end);
    % get GM and WM tissue class
    for ii=1:2
        job_Dnorm.multsdata.vols_tc{ii} = ...
            spm_file(out_US.tiss(ii).c,'number',1);
    end
    % get the parametric maps
    job_Dnorm.multsdata.vols_pm = job.vols_pm;
    % get the warps
    job_Dnorm.multsdata.vols_field = out_Dwarp.files;
    
    % Run the Dartel-Normalize-to-MNI job
    fprintf('\nhMRI-pipeline: running the Darte-normalize-to-MNI module.\n')
    out_Dnorm = hmri_run_proc_dartel_norm(job_Dnorm);
end

% 3/ Deal with smoothing, with hmri_run_proc_smooth
%--------------------------------------------------
proc_smooth = tbx_scfg_hmri_proc_smooth;
[~, job_smooth] = harvest(proc_smooth, proc_smooth, false, false);

% Get the image data, working *only* with mwc1 and mwc2 (GM and WM)
switch job.pipe_c
    case 1 % US+smooth
        job_smooth = fill_fn_from_US(job_smooth,out_US);
    case 2 % US+Dartel+smooth
        % Fit in DARTEL data
        job_smooth.vols_pm = out_Dnorm.vols_wpm;
        job_smooth.vols_tc = out_Dnorm.vols_mwtc;
    otherwise
        error('hmri:pipeline', 'Wrong hmri-pipeline option.');
end
% Get the smoothing kernel
job_smooth.fwhm = job.fwhm;

% Run the *_proc_smooth
fprintf('\nhMRI-pipeline: running the weighted-average (smoothing) module.\n')
out_wa = hmri_run_proc_smooth(job_smooth);
% where the 'out_wa' is organized as a structure out_wa.tc where
% - tc is a cell-array of size {n_TCs x n_pams}
% - each element tc{ii,jj} is a cell array {n_subj x 1} with each subject's
%   smoothed data for the ii^th TC and jj^th MPM


% 4/ Collect output and as needed
out = out_wa; % -> good enouh for the moment!

end

%% ________________________________________________________________________
%
% SUBFUNCTION
%__________________________________________________________________________

function job_smooth = fill_fn_from_US(job_smooth,out_US)
% Fill in the filenames of images (parametric maps and tissue classes) for
% the US+smooth pipeline.

% Parametric maps -> use all of them
N_pm = numel(out_US.maps);
for ii=1:N_pm
    job_smooth.vols_pm{ii} = out_US.maps(ii).wvols_pm;
end
% Tissue classes -> use GM and WM, i.e. #1 and #2
for ii=1:2
    job_smooth.vols_tc{ii} = spm_file(out_US.tiss(ii).mwc,'number',1);
end

% NOTE:
% Not sure it's necessary to add the ',1' to specify the volume for the
% tissue class maps but that's how it looks when using the module with the
% batch GUI -> stay on the safe side and apply.

end
%_______________________________________________________________________
