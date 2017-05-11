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

% % % % % called if DWI processing toolbox not available
% % % % nohmri         = cfg_exbranch;
% % % % nohmri.tag     = 'nohmri';
% % % % nohmri.name    = 'Not available! Help!';
% % % % nohmri.help    = {
% % % %     'No hMRI toolbox available!'   
% % % %     ['The hMRI toolbox does not seem to be available on this computer. ',...   
% % % %     'The directory containing the toolbox implementation should be in the ',...
% % % %     'Matlab path to be used in SPM. Contact Evelyne Balteau if you wish ',...
% % % %     'to use this toolbox (e.balteau@ulg.ac.be).']
% % % %     ['This toolbox should be considered as only a beta (trial) version, ',...
% % % %     'and will include a number of (as yet unspecified) extensions in ',...
% % % %     'future updates. Please report any bugs or problems to e.balteau@ulg.ac.be.']
% % % %     }';
% % % % nohmri.val     = {};
% % % % 
% % % % % redirect to DWI processing toolbox if available, nohmri otherwise
% % % % hmri         = cfg_choice;
% % % % hmri.tag     = 'hMRI';
% % % % hmri.name    = 'hMRI Toolbox';
% % % % hmri.help    = {
% % % %     'This toolbox is a collection of tools necessary for hMRI.'
% % % %     ['This toolbox should be considered as only a beta (trial) version, ',...
% % % %     'and will include a number of (as yet unspecified) extensions in ',...
% % % %     'future updates. Please report any bugs or problems to e.balteau@ulg.ac.be.']
% % % %     }';
% % % % if exist('tbx_cfg_hmri.m','file')
% % % %     hmri.values  = {tbx_cfg_hmri};
% % % % else
% % % %     hmri.values  = {nohmri};
% % % % end
% % % % end
% % % % 
