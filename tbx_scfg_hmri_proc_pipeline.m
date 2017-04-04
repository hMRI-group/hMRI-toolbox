function proc_pipel = tbx_scfg_hmri_proc_pipeline
% Configuration file for the pipeline part of the processing modules of
% the "histological MRI" (hMRI) toolbox.
% -> Provides standard processign pipelines.
% 
% For simplicity, 2 standard pipelines are also set up:
% - US+Smooth -> applies US, warps into MNI, then smoothes
%               (weighted-average)
% - US+Dartel+Smooth -> applies US, builds Dartel template and warps into
%                       MNI, then smoothes (weighted-average)
% Most of the parameters are therefore pre-defined and hardcoded!
% For more flexibility, you ought to use the individual modules and build
% your own pipeline.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

% -------------------------------------------------------------------------
% vols Volumes
% ---------------------------------------------------------------------
vols            = cfg_files;
vols.tag        = 's_vols';
vols.name       = 'Structural images (T1w or MT) for segmentation';
vols.help       = {['Select structural images, i.e. T1w or MT, for ',...
    '"unified segmentation". They are used to create the individuam ',...
    'tissue class maps, e.g. GM and WM posterior probability maps']};
vols.filter     = 'image';
vols.ufilter    = '.*';
vols.num        = [1 Inf];

% ---------------------------------------------------------------------
% vols_pm Parametric maps
% ---------------------------------------------------------------------
vols_pm         = cfg_files;
vols_pm.tag     = 'vols_pm';
vols_pm.name    = 'Parametric maps';
vols_pm.help    = {['Select whole brain parameter maps (e.g. MT, R2*, ',...
    'FA etc) for processing.']};
vols_pm.filter  = 'image';
vols_pm.ufilter = '.*';
vols_pm.num     = [1 Inf];

% ---------------------------------------------------------------------
% Gaussian FWHM
% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Gaussian FWHM';
fwhm.val     = {[6 6 6]};
fwhm.strtype = 'e';
fwhm.num     = [1 3];
fwhm.help    = {['Specify the full-width at half maximum (FWHM) of the ',...
    'Gaussian blurring kernel in mm. Three values should be entered',...
    'denoting the FWHM in the x, y and z directions.']};

% ---------------------------------------------------------------------
% proc_pipel Preprocess maps -> pipelines
% ---------------------------------------------------------------------
proc_pipel         = cfg_exbranch;
proc_pipel.tag     = 'proc_pipel';
proc_pipel.name    = 'Process hMRI maps with pipelines';
proc_pipel.help    = {
    ['Parameter maps are spatially processed and brought into standard space',...
    'for furhter statistical analysis.']
    [' ']
    ['For simplicity, 2 standard pipelines are also set up:']
    ['- US+Smooth -> applies US, warps into MNI, then smoothes (weighted-average)']
    ['US+Dartel+Smooth -> applies US, builds Dartel template and warps' ...
    'into MNI, then smoothes (weighted-average)']
    }'; %#ok<*NBRAK>
proc_pipel.values  = {vols vols_pm fwhm};
proc_pipel.prog = @hmri_run_local_proc_pipeline;
proc_pipel.vout = @vout_preproc;

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% =======================================================================
% VOUT function
% =======================================================================
% TO ADD:
% Need for a check function to ensure the same number of files in each
% series of maps + reference structural.

function dep = vout_preproc(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

cdep = cfg_dep;

% Collect tissue class images (4 of them)
for i=1:numel(job.tissue)
    if job.tissue(i).native(1)
        cdep(end+1) = cfg_dep; %#ok<*AGROW>
        cdep(end).sname = sprintf('c%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'c', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}}); % cfg_findspec({{'filter','nifti'}});
    end
    if job.tissue(i).native(2)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('rc%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'rc', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).warped(1)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('wc%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'wc', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).warped(2)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('mwc%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'mwc', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

% Collect warped parametric maps
for i=1:numel(job.many_sdatas.vols_pm)
    cdep(end+1) = cfg_dep;
    cdep(end).sname = sprintf('Warped par. vols #%d', i);
    cdep(end).src_output = substruct('.', 'maps', '()', {i}, '.', 'wvols_pm', '()', {':'});
    cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%     cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
end

% Collect the deformation fields
cdep(end+1) = cfg_dep;
cdep(end).sname = 'Def. fields';
cdep(end).src_output = substruct('.', 'def', '.', 'fn', '()', {':'});
cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
% cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});

dep = cdep(2:end);
end
