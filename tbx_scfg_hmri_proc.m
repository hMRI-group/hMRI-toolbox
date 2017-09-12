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
    'for furhter statistical analysis.']
    ['This include 3 main processing steps:']
    ['- Unified Segmentation (US) -> produces individual tissue maps' ...
    '+ warping into MNI space']
    ['- DARTEL -> improves the warping into a common space of multiple' ...
    'subjects, based on their tissue maps']
    ['- Smoothing -> tissue specific weighted averaging']
    }'; %#ok<*NBRAK>
proc_modul.values  = {tbx_scfg_hmri_proc_US tbx_scfg_hmri_proc_Dartel tbx_scfg_hmri_proc_smooth};

% ---------------------------------------------------------------------
% proc Preprocess maps
% ---------------------------------------------------------------------
proc         = cfg_choice;
proc.tag     = 'proc';
proc.name    = 'Process hMRI maps';
proc.help    = {
    ['Parameter maps are spatially processed and brought into standard space',...
    'for furhter statistical analysis.']
    ['This include 3 main processing steps:']
    ['- Unified Segmentation (US) -> produces individual tissue maps' ...
    '+ warping into MNI space']
    ['- DARTEL -> improves the warping into a common space of multiple' ...
    'subjects, based on their tissue maps']
    ['- Smoothing -> tissue specific weighted averaging']
    [' ']
    ['For simplicity, 2 standard pipelines are also set up:']
    ['- US+Smooth -> applies US, warps into MNI, then smoothes (weighted-average)']
    ['US+Dartel+Smooth -> applies US, builds Dartel template and warps' ...
    'into MNI, then smoothes (weighted-average)']
    }'; %#ok<*NBRAK>
proc.values  = {tbx_scfg_hmri_proc_pipeline proc_modul};

end
