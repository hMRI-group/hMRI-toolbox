function dicom = tbx_scfg_hmri_dicom_import
% This is a diversion script redirecting towards the DICOM import
% implementation including JSON metadata. The latter is contained in
% 'spm12' and could be used on its own as an extension to SPM12. To do so,
% the 'spm12' directory must be copied into your own SPM12 directory. It
% will then overwrite a few files, but the standard DICOM import
% functionality is preserved whenever the metadata option is disabled.
%
% WARNING: to have the DICOM import working properly from the hMRI-Toolbox
% menu (Batch > SPM > Tools > hMRI-Toolbox > DICOM Import) without
% overwriting SPM12 files, the hMRI scripts must be prioritary in the
% Matlab search path!
%
%_______________________________________________________________________
% Cyclotron Research Centre - University of Liège
% Evelyne Balteau - April 2017
% ======================================================================

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'spm12')));

dicom = spm_cfg_dicom;