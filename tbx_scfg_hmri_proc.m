function proc = tbx_scfg_hmri_proc
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with the processign of the maps
% This include 3 main processing steps:
% - Unified Segmentation (US) -> produces individual tissue maps + warping 
%   into MNI space
% - DARTEL -> improves the warping into a common space of multiple
%   subjects, based on their tissue maps
% - Smoothing -> tissue specific weighted averaging
% These steps are split into 3 sub-modules.
% 
%  For simplicity, 2 standard pipelines are also set up:
% - US+Smooth -> applies US, warps into MNI, then smoothes
%               (weighted-average)
% - US+Dartel+Smooth -> applies US, builds Dartel template and warps into
%                       MNI, then smoothes (weighted-average)
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips

% NOTE:
% Only the Dartel tool is made available. Shoot could be added later on if
% needed or wanted...

% ---------------------------------------------------------------------
% proc_modul Preprocess maps -> individual steps
% ---------------------------------------------------------------------
proc_modul         = cfg_choice;
proc_modul.tag     = 'proc_modul';
proc_modul.name    = 'Proc. hMRI -> Individual modules';
proc_modul.help    = {
    ['Parameter maps are spatially processed and brought into standard space',...
    'for further statistical analysis.']
    ['This includes 4 main processing steps:']
    ['- Unified Segmentation (US) -> produces individual tissue maps' ...
    '+ warping into MNI space']
    ['- DARTEL -> improves the warping into a common space of multiple' ...
    'subjects, based on their tissue maps']
    ['- Smoothing -> tissue specific weighted averaging']
    ['- Mask creation -> tissue specific mask for SPM analysis']
    }'; %#ok<*NBRAK>
proc_modul.values  = {tbx_scfg_hmri_proc_US tbx_scfg_hmri_proc_Dartel ...
    tbx_scfg_hmri_proc_smooth tbx_scfg_hmri_proc_crtMask};

% ---------------------------------------------------------------------
% proc Preprocess maps
% ---------------------------------------------------------------------
proc         = cfg_choice;
proc.tag     = 'proc';
proc.name    = 'Process hMRI maps';
proc.help    = {
    ['Parameter maps are spatially processed and brought into standard space',...
    'for further statistical analysis.']
    ['This includes 4 main processing steps:']
    ['- Unified Segmentation (US) -> produces individual tissue maps' ...
    '+ standard warping into MNI space']
    ['- DARTEL -> improves the warping into a population common space, ' ...
    'based on the subjects'' GM and WM tissue maps']
    ['- Smoothing -> tissue specific weighted averaging']
    ['- Mask creation -> creating tissue specific masks']
    [' ']
    ['For simplicity, 2 standard pipelines are also set up:']
    ['- US+Smooth+MaskCrt -> applies US, warps into MNI, smoothes (weighted-average), then creates masks']
    ['- US+Dartel+Smooth+MaskCrt -> applies US, builds Dartel template and warps' ...
    'into MNI, smoothes (weighted-average), then creates masks']
    }'; %#ok<*NBRAK>
proc.values  = {tbx_scfg_hmri_proc_pipeline proc_modul};

end
