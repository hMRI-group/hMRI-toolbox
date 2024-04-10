function [denoisingtype,DNparameters] = tbx_scfg_hmri_denoising_menu

%--------------------------------------------------------------------------
% Denoising (DN) defaults file
%--------------------------------------------------------------------------
DNdefaults         = cfg_files;
DNdefaults.tag     = 'DNdefaults';
DNdefaults.name    = 'Customised DN defaults file';
DNdefaults.help    = {['Select the [hmri_denoising*_defaults_*.m] file containing ' ...
    'the parameters to process the denoising.'], ...
    ['Please make sure that the modified denoising (optional) parameters ' ...
    'are correct for your data. To create your own customised defaults file, ' ...
    'edit the distributed version and save it with a meaningful name such as ' ...
    '[hmri_denoising_local_defaults.m].']};
DNdefaults.filter  = 'm';
DNdefaults.dir     = fullfile(fileparts(mfilename('fullpath')),'config','local');
DNdefaults.ufilter = [...
    '(?i)',   ... ignore case
    '^hmri_', ... filename starts with 'hmri_'
    '.*denoising.*', ... filename includes 'denoising' somewhere in the middle
    '\.m$'    ... filename ends with '.m'
    ];
DNdefaults.num     = [1 1];

% ---------------------------------------------------------------------
% Use metadata or standard defaults (no customization)
% ---------------------------------------------------------------------
DNmetadata         = cfg_entry;
DNmetadata.tag     = 'DNmetadata';
DNmetadata.name    = 'Use metadata or standard defaults';
DNmetadata.help    = {''};
DNmetadata.strtype = 's';
DNmetadata.num     = [1 Inf];
DNmetadata.val     = {'yes'};

%--------------------------------------------------------------------------
% DN processing parameters
%--------------------------------------------------------------------------
DNparameters        = cfg_choice;
DNparameters.tag    = 'DNparameters';
DNparameters.name   = 'Processing parameters';
DNparameters.help   = {['You can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own customised defaults file. ' ...
    'For possible modification of default parameters, please be careful about which parameters are '...
    'required and which parameters are optional.']};
DNparameters.values = {DNmetadata DNdefaults};
DNparameters.val    = {DNmetadata};

% ---------------------------------------------------------------------
% Magnitude input images 
% ---------------------------------------------------------------------
mag_img         = cfg_files;
mag_img.tag     = 'mag_input';
mag_img.name    = 'Magnitude input (required)';
mag_img.help    = {'Select the (required) magnitude images to be denoised'};
mag_img.filter  = 'image';
mag_img.ufilter = '.*';
mag_img.num     = [1 Inf];


% ---------------------------------------------------------------------
% Phase input images
% ---------------------------------------------------------------------
phase_img         = cfg_files;
phase_img.tag     = 'phase_input';
phase_img.name    = 'Phase input (optional)';
phase_img.help    = {'Select the (optional) phase images to be denoised'};
phase_img.filter  = 'image';
phase_img.ufilter = '.*';
phase_img.num     = [0 Inf];

% ---------------------------------------------------------------------
% Standard deviation parameter
% ---------------------------------------------------------------------
std         = cfg_entry;
std.tag     = 'std';
std.name    = 'Standard deviation cut-off';
std.val     = {[1.05]};
std.strtype = 'e';
std.num     = [1 1];
std.help    = {['Specify the standard deviation cut off for denoising']};


% ---------------------------------------------------------------------
% Neighborhood size
% ---------------------------------------------------------------------
ngbsize         = cfg_entry;
ngbsize.tag     = 'ngbsize';
ngbsize.name    = 'Neighborhood size';
ngbsize.val     = {[4]};
ngbsize.strtype = 'e';
ngbsize.num     = [1 1];
ngbsize.help    = {['Specify the neghborhood size']};

% ---------------------------------------------------------------------
% LCPCA Denoising protocol
% ---------------------------------------------------------------------
denoisinginput_lcpca      = cfg_branch;
denoisinginput_lcpca.tag  = 'lcpca_denoise';
denoisinginput_lcpca.name = 'LCPCA denoising';
denoisinginput_lcpca.help = {'Input Magnitude/Phase images for Lcpca-denoising'
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_denoisinglocal_defaults_*.m] customised defaults file.']};
denoisinginput_lcpca.val  = {mag_img phase_img DNparameters std ngbsize};


% ---------------------------------------------------------------------
% menu denoisingtype
% ---------------------------------------------------------------------
denoisingtype        = cfg_choice;
denoisingtype.tag    = 'denoisingtype';
denoisingtype.name   = 'Denoising method for raw/processed images';
denoisingtype.help   = {'Choose the methods for denoising.'
    ['Various types of denoising methods can be handled by the hMRI ' ...
    'toolbox. See list below for a ' ...
    'brief description of each type. Note that all types may not be ' ...
    'available at your site.']
    ['- Lcpca denoising: Bazin, et al. (2019) Denoising High-Field Multi-Dimensional MRI With Local'...
'Complex PCA, Front. Neurosci. doi:10.3389/fnins.2019.01066']
    [' - No denoising:...']};
denoisingtype.values = {denoisinginput_lcpca};
denoisingtype.val    = {denoisinginput_lcpca};


end