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
% vols_tc Parametric volumes
% ---------------------------------------------------------------------
vols_tc         = cfg_files;
vols_tc.tag     = 'vols_tc';
vols_tc.name    = 'mwTC images';
vols_tc.help    = {'Select the modulated warped tissue classes (mwc*)'};
vols_tc.filter  = 'image';
vols_tc.ufilter = '^mwc.*';
vols_tc.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_TCs Tissue class (TC) maps, used for 'many subjects'
% ---------------------------------------------------------------------
m_TCs            = cfg_repeat;
m_TCs.tag        = 'maps';
m_TCs.name       = 'Modulated warped tissue class';
m_TCs.values     = {vols_tc };
m_TCs.val        = {vols_tc };
m_TCs.num = [1 Inf];
m_TCs.help       = {['Select the modulated warped tissue classes (TC) ',...
    'of interest from all subjects. This ould typically be the mwc1* ',...
    'and mwc2* images for GM and WM.']};

% ---------------------------------------------------------------------
% tpm Tissue Probability Maps
% ---------------------------------------------------------------------
tpm         = cfg_files;
tpm.tag     = 'tpm';
tpm.name    = 'Tissue probability maps';
tpm.help    = {'Select the TPM used for the segmentation.'};
tpm.filter  = 'image';
tpm.ufilter = '.*';
tpm.num     = [1 1];
tpm.val     = {{fullfile(spm('dir'),'toolbox','hMRI','tpm','unwTPM_sl2.nii,1')}};

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

%% EXEC function
% ---------------------------------------------------------------------
% proc_smooth Processing hMRI -> smoothing
% ---------------------------------------------------------------------
proc_smooth         = cfg_exbranch;
proc_smooth.tag     = 'proc_smooth';
proc_smooth.name    = 'Proc. hMRI -> Smoothing';
proc_smooth.val     = {m_pams m_TCs tpm fwhm};
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

n_pams = numel(job.vols_pm);     % #parametric image types
n_TCs = numel(job.vols_tc);      % #tissue classes
% n_subj = numel(job.vols_pm{1});  % #subjects

cdep = cfg_dep;
for ii=1:n_TCs
    for jj=1:n_pams
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('TC #%d, pMap #%d', ii, jj);
        cdep(end).src_output = substruct('.', 'tc', '()', {ii}, '.', 'map', '()', {jj}, '.', 'fn');
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end
    
dep = cdep(2:end);

end

% Check number of files match
function chk = check_proc_smooth(job)
% ensure they are the same for each list of files, one of each per subject.

n_pams = numel(job.vols_pm);
n_TCs = numel(job.vols_tc);
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
    ni_TC = numel(job.vols_tc{1});
    for ii=2:n_TCs
        if ni_TC ~= numel(job.vols_tc{ii})
            chk = [chk 'Incompatible number of TCs. ']; %#ok<*AGROW>
            break
        end
    end
end
if n_pams>0 && n_TCs>0
    if numel(job.vols_pm{1}) ~= numel(job.vols_tc{1})
        chk = [chk 'Incompatible number of maps & TCs.'];
    end
end

end

