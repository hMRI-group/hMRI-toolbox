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

% ---------------------------------------------------------------------
% indir Input directory as output directory
% ---------------------------------------------------------------------
indir         = cfg_menu;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {['Output files will be written to the same folder as ',...
    'each corresponding input file.']};
indir.labels  = {'Yes'};
indir.values  = {1};
indir.val     = {1};

% ---------------------------------------------------------------------
% outdir Output directory for all data
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory, all together';
outdir.help    = {['Select a directory where all output files from all '... 
    'subjects put together will be written to.']};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

% ---------------------------------------------------------------------
% outdir_ps Output directory for per-subject organisation
% ---------------------------------------------------------------------
outdir_ps         = cfg_files;
outdir_ps.tag     = 'outdir_ps';
outdir_ps.name    = 'Output directory, with per-subject sub-directory';
outdir_ps.help    = {['Select a directory where output files will be '...
    'written to, in each subject''s sub-directory.']};
outdir_ps.filter = 'dir';
outdir_ps.ufilter = '.*';
outdir_ps.num     = [1 1];

% ---------------------------------------------------------------------
% output Output choice
% ---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {['Output directory can be the same as the input ',...
    'directory for each input file or user selected (one for everything ',...
    'or preserve a per-subject organisation).']};
output.values  = {indir outdir outdir_ps };
output.val     = {indir};

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
vols_pm.name    = 'Maps (single type)';
vols_pm.help    = {['Select whole brain parameter maps (e.g. MT, R2*, ',...
    'FA, etc.) from all subjects for processing.']};
vols_pm.filter  = 'image';
vols_pm.ufilter = '.*';
vols_pm.num     = [1 Inf];

% ---------------------------------------------------------------------
% many_pams Parameter maps
% ---------------------------------------------------------------------
% used for 'many subjects', i.e. list the data per map type across subjects
many_pams            = cfg_repeat;
many_pams.tag        = 'maps';
many_pams.name       = 'Parametric maps';
many_pams.values     = {vols_pm };
many_pams.val        = {vols_pm };
many_pams.num = [1 Inf];
many_pams.help       = {['Select whole brain parameter maps (e.g. MT, ',...
    'R2*, FA, etc.) from all subjects for processing, one type at a time.']};

% ---------------------------------------------------------------------
% pipe_c Pipeline choice
% ---------------------------------------------------------------------
pipe_c        = cfg_menu;
pipe_c.tag    = 'pipe_c';
pipe_c.name   = 'Pipeline';
pipe_c.help   = {
    'Chose the predefined pipeline that you prefer:'
    '- US+Smooth -> applies US, warps into MNI, then smoothes (weighted-average)'
    ['- US+Dartel+Smooth -> applies US, builds Dartel template and warps into' ...
      'MNI, then smoothes (weighted-average)']
    }';
pipe_c.labels = {
                 'US+smooth'
                 'US+Dartel+smooth'}';
pipe_c.values = {1 2};
pipe_c.val    = {1};

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
proc_pipel.name    = 'Proc. hMRI -> Pipelines';
proc_pipel.help    = {
    ['Parameter maps are spatially processed and brought into standard space',...
    'for furhter statistical analysis.']
    [' ']
    ['For simplicity, 2 standard pipelines are also set up:']
    ['- US+Smooth -> applies US, warps into MNI, then smoothes (weighted-average)']
    ['US+Dartel+Smooth -> applies US, builds Dartel template and warps' ...
    'into MNI, then smoothes (weighted-average)']
    }'; %#ok<*NBRAK>
proc_pipel.val  = {output vols many_pams fwhm pipe_c};
proc_pipel.prog = @hmri_run_proc_pipeline;
proc_pipel.vout = @vout_proc_pipeline;
proc_pipel.check = @check_data;

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% =======================================================================
% VOUT function
% =======================================================================

% Collect and prepare output
function dep = vout_proc_pipeline(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.
% There should be one series of images per parametric map and tissue class,
% e.g. in the usual case of 4 MPMs and GM/WM -> 8 series of image

n_pams = numel(job.vols_pm); % #parametric image types
n_TCs  = 2;                  % #tissue classes = 2, by default

cdep = cfg_dep;
for ii=1:n_TCs
    for jj=1:n_pams
        cdep(end+1) = cfg_dep; %#ok<*AGROW>
        cdep(end).sname = sprintf('TC #%d, pMap #%d', ii, jj);
        cdep(end).src_output = substruct('.', 'tc', '{}', {ii,jj});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end
    
dep = cdep(2:end);

end

%% =======================================================================
% CHECKING the data
% ========================================================================
function t = check_data(job)
% Checking that the data are consistent.
t   = {};

nSubj = numel(job.s_vols); % number of subjects from #struct images
nPara = numel(job.vols_pm); % number of maps type
% Check number of structurals matches the number of parametric maps per
% type
if nPara>0
    for ii=1:nPara
        if numel(job.vols_pm{ii})~=0 && numel(job.vols_pm{ii})~=nSubj
        t{1} = 'Number of maps not matching number of structural images/subjects!';
            warndlg(t,'Maps numbers');
            return
        end    
    end
end

end
