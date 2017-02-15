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

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','hmri')); end

% Work is split into 2 main branches:
% - creating the MPM's from the raw images -> tbx_cfg_hmri_crm
% - spatially processing the MPM's -> tbx_cfg_hmri_proc

% ---------------------------------------------------------------------
% hmri hMRI Tools
% ---------------------------------------------------------------------
hmri         = cfg_choice;
hmri.tag     = 'hmri';
hmri.name    = 'hMRI Tools';
hmri.help    = {
    ['This toolbox is based around the ``Regional specificity of MRI ',...
    'contrast ... (VBQ)'' paper by Draganski et al., 2011 NeuroImage ',...
    'and ''Unified segmentation based correction... (UNICORT) paper by ',...
    'Weiskopf et al., 2011. ']
    ['This toolbox should be considered as only a beta (trial) version, ',...
    'and will include a number of (as yet unspecified) extensions in ',...
    'future updates.  Please report any bugs or problems to lreninfo@gmail.com.']
    }';
hmri.values  = {tbx_scfg_hmri_crm tbx_scfg_hmri_proc };
end

