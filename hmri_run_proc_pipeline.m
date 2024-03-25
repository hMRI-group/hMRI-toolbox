function out = hmri_run_proc_pipeline(job)
% Deal with the preprocessing pipelines. There are 2 options: 
% 1/ "US+Smooth+MaskCrt"
% 2/ "US+Dartel+Smooth+MaskCrt" 
% The only difference is the inclusion of Dartel to bring all the data into 
% a common space, before warping into MNI.
% 
% So the operations go through 3 or 4 main steps:
% 1/ unified segementation, if no Dartel with warping to MNI
% 2/ if with Dartel, Dartel estimation and warping to MNI
% 3/ smoothing with tissue class specific smoothing
% 4/ creation of tissue specific mask
%
% Input include only some parametric maps, the structural maps (for
% segmentation), the required smoothing and which pipeline to use. All
% other options are hard-coded!
% By default, these pipelines focus only on the first 2 tissue classes,
% i.e. GM and WM only.
%
% For more flexibility, individual modules can be combined. :-)
%
% The output structure 'out' provides:
% .tc     : cell-array of size {n_TCs x n_pams}. Each element tc{ii,jj} is 
%           a cell array {n_subj x 1} with each subject's smoothed data for
%           the ii^th TC and jj^th MPM
% .smwc   : cell-array of size {n_TCs x1}. Each element smwc{ii} is a 
%           char array (n_subj x 1) with each subject's smooth modulated 
%           warped ii^th tissue class
% .maskTC : cell-array of size {n_TCs x1}. Each element smwc{ii} is the 
%           filename of the ii^th tissue specific masks
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

%% 1/ Setup the smoothing and execute hmri_run_proc_US
%----------------------------------------------------
% Get the proc_US job structure, with all the defaults
proc_us = tbx_scfg_hmri_proc_US;
[~, job_US] = harvest(proc_us, proc_us, false, false);

% Fill in the data now: parametric maps & structurals for segmentation
if job.pipe_c==1 % US+smooth -> warp the parametric maps here.
    job_US.many_sdatas.vols_pm = job.vols_pm;
end
job_US.many_sdatas.channel.vols = job.s_vols;
% Set voxel sizes & BB
job_US.many_sdatas.vox = job.vox;
job_US.many_sdatas.bb  = job.bb;

% Get the output direcotry across
job_US.many_sdatas.output = job.output;

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

%% 2/ Proceed with dartel (only if requested)
%-------------------------------------------
% including template create and warping into MNI space
if job.pipe_c == 2
    [~,ver] = spm('Ver');
    if str2double(ver)>=7219
        % Knows how to handle output specification
        use_spm_output_handling = true;
    else
        use_spm_output_handling = false;
    end
    
    % DARTEL processing:
    proc_Dartel = tbx_scfg_hmri_proc_Dartel;
    % a) warping with template creation
    proc_Dwarp = proc_Dartel.values{1};
    [~, job_Dwarp] = harvest(proc_Dwarp, proc_Dwarp, false, false);
    
    % Fill in the filenames for rc1 and rc2, i.e. *only* GM and WM
    for ii=1:2
        job_Dwarp.images{ii} = spm_file(out_US.tiss(ii).rc,'number',1);
    end
    job_Dwarp.output = job.output;
    % Run the Dartel-warp job
    fprintf('\nhMRI-pipeline: running the Dartel-warp module.\n')
    out_Dwarp = hmri_run_proc_dartel_template(job_Dwarp);
    %     out_Dwarp = spm_dartel_template(job_Dwarp);
    
    if ~use_spm_output_handling
        % Move if using specific per-subject subdirectory -> 'dart_files'
        % This should also include the flow fields!
        if isfield(job.output,'outdir_ps')
            dn_dartel = fullfile(job.output.outdir_ps{1},'Dartel_Templates');
            if ~exist(dn_dartel,'dir')
                mkdir(dn_dartel)
            end
            current_path = spm_file(out_Dwarp.template{1},'path');
            f2move = spm_select('FPList',current_path,'^Template_[\d]\.nii$');
            for ii=1:size(f2move,1)
                movefile(deblank(f2move(ii,:)),dn_dartel);
            end
            out_Dwarp.template = spm_file(out_Dwarp.template,'path',dn_dartel);
        end
    end
    
    % b) normalize to  MNI
    proc_Dnorm = proc_Dartel.values{3};
    [~, job_Dnorm] = harvest(proc_Dnorm, proc_Dnorm, false, false);
    % get last tempalte
    job_Dnorm.template = out_Dwarp.template(end);
    % get GM, WM and CSF tissue class
    for ii=1:3
        job_Dnorm.multsdata.vols_tc{ii} = ...
            spm_file(out_US.tiss(ii).c,'number',1);
    end
    % get the parametric maps
    job_Dnorm.multsdata.vols_pm = job.vols_pm;
    % get the warps
    job_Dnorm.multsdata.vols_field = out_Dwarp.files;
    % Set voxel sizes
    job_Dnorm.vox = job.vox;
    job_Dnorm.bb  = job.bb;
    % get the output directory
    job_Dnorm.output = job.output;
    
    % Run the Dartel-Normalize-to-MNI job
    fprintf('\nhMRI-pipeline: running the Darte-normalize-to-MNI module.\n')
    out_Dnorm = hmri_run_proc_dartel_norm(job_Dnorm);
end

%% 3/ Deal with smoothing, with hmri_run_proc_smooth
%--------------------------------------------------
proc_smooth = tbx_scfg_hmri_proc_smooth;
[~, job_smooth] = harvest(proc_smooth, proc_smooth, false, false);

% Get the image data, working with mwc1, mwc2 and mwc3 (GM, WM and CSF)
switch job.pipe_c
    case 1 % US+smooth
        job_smooth = fill_fn_from_US(job_smooth,out_US);
    case 2 % US+Dartel+smooth
        % Fit in DARTEL data
        job_smooth.vols_pm = out_Dnorm.vols_wpm;
        job_smooth.vols_mwc = out_Dnorm.vols_mwc;
    otherwise
        error('hmri:pipeline', 'Wrong hmri-pipeline option.');
end
% Get the smoothing kernel
job_smooth.fwhm = job.fwhm;

% Run the *_proc_smooth
fprintf('\nhMRI-pipeline: running the weighted-average (smoothing) module.\n')
out_wa = hmri_run_proc_smooth(job_smooth);
% The 'out_wa' structure is organized as a structure with 2 fields
% .tc   : cell-array of size {n_TCs x n_pams}. Each element tc{ii,jj} is a 
%         cell array {n_subj x 1} with each subject's smoothed data for
%         the ii^th TC and jj^th MPM
% .smwc : cell-array of size {n_TCs x1}. Each element smwc{ii} is a char
%         array (n_subj x 1) with each subject's smooth modulated warped
%         ii^th tissue class

%% 4/ Create the tissue specific masks
proc_crtMask = tbx_scfg_hmri_proc_crtMask;
[~, job_crtMask] = harvest(proc_crtMask, proc_crtMask, false, false);
% Get the smwc images in & fix output directory
job_crtMask.vols_smwc = out_wa.smwc;
if isfield(job.output,'outdir_ps')
    job_crtMask.options.output = struct('outdir',job.output.outdir_ps);
else
    job_crtMask.options.output = job.output;
end
% Call to 
out_cM = hmri_run_proc_crtMask(job_crtMask);
% The out_cM structure has 2 fields:
% .fn_maskTC : filenames (char array) of the tissue specific masks
% .fn_meanTC : filenames (char array) of the mean tissue class images


%% 5/ Collect output and as needed
out = out_wa; 
% There are 2 fields with cell array in 'out_wa' :
% .tc   : tissue specifically smoothed qMRIs 
% .smwc : smooth modulated warped tissue class images

% Adding the the tissue mask images -> necessary for SPM analysis
out.maskTC = cellstr(out_cM.fn_maskTC);

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
% Tissue classes -> use GM, WM and CSF, i.e. #1, #2 and #3
for ii=1:3
    job_smooth.vols_mwc{ii} = spm_file(out_US.tiss(ii).mwc,'number',1);
end

% NOTE:
% Not sure it's necessary to add the ',1' to specify the volume for the
% tissue class maps but that's how it looks when using the module with the
% batch GUI -> stay on the safe side and apply.

end
%_______________________________________________________________________
