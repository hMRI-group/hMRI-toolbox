function hmri = tbx_cfg_hmri_redirect
% Diversion config file for the hMRI toolbox
%
% PURPOSE
% The present hMRI config file redirects towards the full implementation of
% the toolbox (tbx_cfg_hmri) only if present in the Matlab path. The
% toolbox can therefore be stored in a directory independent from the SPM
% implementation and synch'd with the main hMRI-Toolbox repository whenever
% needed. If tbx_cfg_hmri is not found in the Matlab path, the hMRI Tools
% are listed in the SPM Batch GUI but not available. A brief help section
% provides the user with instructions for hMRI Tools installation.
%
% USAGE
% Copy this file into a directory in the SPM toolbox directory (e.g.
% <path-to-SPM>/toolbox/hMRI-Redirect). Add the hMRI-Toolbox directory
% (containing the full implementation of the toolbox) to the Matlab path.
% Restart SPM and the Batch GUI. The hMRI Tools will be available in the
% SPM>Tools menu.
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%__________________________________________________________________________
% Cyclotron Research Centre - University of Liège
% Evelyne Balteau - April 2017
%==========================================================================

if ~isdeployed, addpath(fileparts(mfilename('fullpath'))); end

try
    % hMRI is available
    hmri = tbx_cfg_hmri;
catch %#ok<CTCH>
    % No hMRI toolbox found
    hmri         = cfg_exbranch;
    hmri.tag     = 'hMRI';
    hmri.name    = 'hMRI Tools - not available!';
    hmri.help    = {
        ['The hMRI toolbox does not seem to be available on this computer. ',...
        'The directory containing the toolbox implementation should be ',...
        'in the Matlab path to be used in SPM. See installation ',...
        'instructions on the hMRI-Toolbox repository: ']
        'https://github.molgen.mpg.de/VBQ-toolbox-group/hMRI-Toolbox/wiki.'
        }';
    hmri.val  = {};
end

end

