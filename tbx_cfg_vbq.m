function vbq = tbx_cfg_vbq
% Configuration file for the Voxel-Based Quantification (VBQ)
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%_______________________________________________________________________
%
% Bogdan Draganski & Ferath Kherif, 2011
% ======================================================================

% $Id$

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','VBQ')); end

% Work is split into 2 main branches:
% - creating the MPM's from the raw images -> tbx_cfg_vbq_crm
% - spatially processing the MPM's -> tbx_cfg_vbq_proc

% ---------------------------------------------------------------------
% vbq VBQ Tools
% ---------------------------------------------------------------------
vbq         = cfg_choice;
vbq.tag     = 'VBQ';
vbq.name    = 'VBQ Tools';
vbq.help    = {
    ['This toolbox is based around the ``Regional specificity of MRI ',...
    'contrast ... (VBQ)'' paper by Draganski et al., 2011 NeuroImage ',...
    'and ''Unified segmentation based correction... (UNICORT) paper by ',...
    'Weiskopf et al., 2011. ']
    ['This toolbox should be considered as only a beta (trial) version, ',...
    'and will include a number of (as yet unspecified) extensions in ',...
    'future updates.  Please report any bugs or problems to lreninfo@gmail.com.']
    }';
vbq.values  = {tbx_cfg_vbq_crm tbx_cfg_vbq_proc };
end

