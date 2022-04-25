function hmri = tbx_cfg_hmri
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%_______________________________________________________________________
%
% Bogdan Draganski & Ferath Kherif, 2011
% ======================================================================

if ~isdeployed, 
    hMRIpath = fileparts(mfilename('fullpath'));
    addpath(hMRIpath);
    addpath(fullfile(hMRIpath, 'config'));
    addpath(fullfile(hMRIpath, 'etpm'));
    addpath(fullfile(hMRIpath, 'spm12'));
    addpath(fullfile(hMRIpath, 'spm12','config'));
    addpath(fullfile(hMRIpath, 'spm12','metadata'));
    addpath(genpath(fullfile(hMRIpath, 'EPG_for_hmri_tbx')));
end

% The toolbox is split into 5 main modules:
% - Configure toolbox -> tbx_scfg_hmri_config
% - DICOM import -> tbx_scfg_hmri_dicom_import
% - Auto-reorient -> tbx_scfg_hmri_autoreorient
% - Create hMRI maps -> tbx_scfg_hmri_create
% - Process (spatially) hMRI maps -> tbx_scfg_hmri_proc

% ---------------------------------------------------------------------
% hmri hMRI Tools
% ---------------------------------------------------------------------
hmri         = cfg_choice;
hmri.tag     = 'hmri';
hmri.name    = 'hMRI Tools';
hmri.help    = {
    ['This toolbox is based around the ''Regional specificity of MRI ',...
    'contrast ... (VBQ)'' paper by Draganski et al., 2011 NeuroImage ',...
    'and ''Unified segmentation based correction... (UNICORT)'' paper by ',...
    'Weiskopf et al., 2011. ']
    ['This toolbox should be considered as only a beta (trial) version, ',...
    'and will include a number of (as yet unspecified) extensions in ',...
    'future updates.  Please report any bugs or problems to lreninfo@gmail.com.']
    }';
hmri.values  = {tbx_scfg_hmri_config tbx_scfg_hmri_dicom_import tbx_scfg_hmri_autoreorient tbx_scfg_hmri_imperf_spoil tbx_scfg_hmri_B1_create tbx_scfg_hmri_create tbx_scfg_hmri_wcomb tbx_scfg_hmri_quality tbx_scfg_hmri_proc tbx_scfg_hmri_QUIQI};
end

