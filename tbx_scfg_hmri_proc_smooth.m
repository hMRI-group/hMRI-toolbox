function proc_smooth = tbx_scfg_hmri_proc_smooth
% Configuration file for the "smoothing", weighted averaging, of
% quantitative maps
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Center

% Written by Christophe Phillips

% NOTE:
% data are selected in a 'many subject' style, i.e. all the images of one
% type are selected from many subjects at once!
% 
% It could be advantageous to define the TPM in a definition file and use
% it when ever we need it. Right now, this is hard-coded in the cfg file!
% Same goes for a bunch of parameters for each tissue class (e.g. number of
% Guassians, what is written out, bias correction, etc.) which ar enow +/-
% hard-coded in this file.

% ---------------------------------------------------------------------
% vols_pm Parametric volumes
% ---------------------------------------------------------------------
vols_pm         = cfg_files;
vols_pm.tag     = 'vols_pm';
vols_pm.name    = 'Volumes';
vols_pm.help    = {['Select whole brain parameter maps (e.g. MT, R2*, ',...
    'FA etc) warped into MNI space.']};
vols_pm.filter  = 'image';
vols_pm.ufilter = '^w.*';
vols_pm.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_pams Parameter maps, used for 'many subjects'
% ---------------------------------------------------------------------
m_pams            = cfg_repeat;
m_pams.tag        = 'm_pams';
m_pams.name       = 'Warped parameter maps';
m_pams.values     = {vols_pm };
m_pams.val        = {vols_pm };
m_pams.num = [1 Inf];
m_pams.help       = {['Select whole brain parameter maps (e.g. MT, ',...
    'R2*, FA etc) warped into MNI space.']};

% ---------------------------------------------------------------------
% vols_mwc Modulated warped tissue segement volumes
% ---------------------------------------------------------------------
vols_mwc         = cfg_files;
vols_mwc.tag     = 'vols_mwc';
vols_mwc.name    = 'mwc images';
vols_mwc.help    = {'Select the modulated warped tissue segements (mwc*).', ... 
    'Pick only one type of mwc* images across all subjects!.'};
vols_mwc.filter  = 'image';
vols_mwc.ufilter = '^mwc.*';
vols_mwc.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_MWCs Modulate warped tissue segement (MWC) maps
% ---------------------------------------------------------------------
m_MWCs            = cfg_repeat;
m_MWCs.tag        = 'm_MWCs';
m_MWCs.name       = 'Modulated warped tissue segements';
m_MWCs.values     = {vols_mwc };
m_MWCs.val        = {vols_mwc };
m_MWCs.num = [1 Inf];
m_MWCs.help       = {['Select the modulated warped tissue segments ',...
    'of interest from all subjects.'], ...
    ['For the typical case of GM and WM, you would selectall the mwc1* images ', ...
    'in one set of ''mwc_images'' and the mwc2* ones in second set of ', ...
    '''mwc_images''!']};

% ---------------------------------------------------------------------
% tpm Tissue Probability Maps
% ---------------------------------------------------------------------
% use the hMRI specific TPMs.
fn_tpm = hmri_get_defaults('proc.TPM');

tpm         = cfg_files;
tpm.tag     = 'tpm';
tpm.name    = 'Tissue probability maps';
tpm.help    = {'Select the TPM used for the segmentation.'};
tpm.filter  = 'image';
tpm.ufilter = '.*';
tpm.num     = [1 1];
tpm.val     = {{fn_tpm}};

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

%% EXEC function
% ---------------------------------------------------------------------
% proc_smooth Processing hMRI -> smoothing
% ---------------------------------------------------------------------
proc_smooth         = cfg_exbranch;
proc_smooth.tag     = 'proc_smooth';
proc_smooth.name    = 'Proc. hMRI -> Smoothing';
proc_smooth.val     = {m_pams m_MWCs tpm fwhm output};
proc_smooth.check   = @check_proc_smooth;
proc_smooth.help    = { 
    'Applying tissue specific smoothing, aka. weighted averaging, ', ...
    'in order to limit partial volume effect.'};
proc_smooth.prog = @hmri_run_proc_smooth;
proc_smooth.vout = @vout_smooth;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect and prepare output
function dep = vout_smooth(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.
% There should be one series of images per parametric map and tissue class,
% e.g. in the usual case of 4 MPMs and GM/WM -> 8 series of image
% + 1 series of smoothed tissue class images

n_pams = numel(job.vols_pm);     % #parametric image types
n_TCs = numel(job.vols_mwc);     % #tissue classes

cdep = cfg_dep;
for ii=1:n_TCs
    for jj=1:n_pams
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('TC #%d, pMap #%d', ii, jj);
        cdep(end).src_output = substruct('.', 'tc', '{}', {ii,jj});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end
    
for ii=1:n_TCs
    cdep(end+1) = cfg_dep;
    cdep(end).sname = sprintf('smwTC #%d', ii);
    cdep(end).src_output = substruct('.', 'smwc', '{}', {ii});
    cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
end
dep = cdep(2:end);

end

% Check number of files match
function chk = check_proc_smooth(job)
% ensure they are the same for each list of files, one of each per subject.

n_pams = numel(job.vols_pm);
n_TCs = numel(job.vols_mwc);
chk = '';

if n_pams>1
    ni_pams = numel(job.vols_pm{1});
    for ii=2:n_pams
        if ni_pams ~= numel(job.vols_pm{ii})
            chk = [chk 'Incompatible number of maps. '];
            break
        end
    end
end
if n_TCs>1
    ni_TC = numel(job.vols_mwc{1});
    for ii=2:n_TCs
        if ni_TC ~= numel(job.vols_mwc{ii})
            chk = [chk 'Incompatible number of tissue segments. ']; %#ok<*AGROW>
            break
        end
    end
end
if n_pams>0 && n_TCs>0
    if numel(job.vols_pm{1}) ~= numel(job.vols_mwc{1})
        chk = [chk 'Incompatible number of maps & tissue segments.'];
    end
end

end

