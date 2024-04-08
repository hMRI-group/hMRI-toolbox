function [denoisingtype,DNparameters] = tbx_scfg_hmri_denoising_menu

%--------------------------------------------------------------------------
% Denoising (DN) defaults file
%--------------------------------------------------------------------------
DNdefaults         = cfg_files;
DNdefaults.tag     = 'DNdefaults';
DNdefaults.name    = 'Customised DN defaults file';
DNdefaults.help    = {['Select the [hmri_denoising*_defaults_*.m] file containing ' ...
    'the parameters to process the denoising By default, parameters will be ' ...
    'collected from metadata when available. Defaults parameters are provided as ' ...
    'fallback solution when metadata are not available and/or uncomplete.'], ...
    ['Please make sure that the parameters defined in the defaults file ' ...
    'are correct for your data. To create your own customised defaults file, ' ...
    'edit the distributed version and save it with a meaningful name such as ' ...
    '[hmri_denoising_local_defaults_*myprotocol*.m].']};
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
    'defaults parameters (recommended) or select your own customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
DNparameters.values = {DNmetadata DNdefaults};
DNparameters.val    = {DNmetadata};
end