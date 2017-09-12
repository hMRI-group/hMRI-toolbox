function hmri_config = tbx_scfg_hmri_config
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with local defaults ("Configure Toolbox" branch)
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging


%--------------------------------------------------------------------------
% Customised defaults
%--------------------------------------------------------------------------
customised         = cfg_files;
customised.tag     = 'customised';
customised.name    = 'Customised';
customised.help    = {['Select the [hmri_local_defaults_*.m] file containing ' ...
    'the specific defaults to process your data. Note that all other defaults ' ...
    'values will be reinitialised to their standard values.']};
customised.filter  = 'm';
customised.dir     = fullfile(fileparts(mfilename('fullpath')),'config','local');
customised.ufilter = '^hmri_.*\.m$';
customised.num     = [1 1];
%customised.def     = @(val)hmri_get_defaults('local_defaults', val{:});

% ---------------------------------------------------------------------
% Use standard defaults (no customization)
% ---------------------------------------------------------------------
standard           = cfg_entry;
standard.tag       = 'standard';
standard.name      = 'Standard';
standard.help      = {''};
standard.strtype = 's';
standard.num     = [1 Inf];
standard.val     = {'yes'};

%--------------------------------------------------------------------------
% Set hMRI defaults parameters
%--------------------------------------------------------------------------
hmri_setdef         = cfg_choice;
hmri_setdef.tag     = 'hmri_setdef';
hmri_setdef.name    = 'Defaults parameters';
hmri_setdef.help    = {['You can either stick with standard defaults parameters ' ...
    'from [hmri_defaults.m] or select your own customised defaults file.']};
hmri_setdef.values  = {standard customised};
hmri_setdef.val     = {standard};


% ---------------------------------------------------------------------
% Configure the hMRI Toolbox - load local, user-defined defaults file and
% overwrite standard defaults 
% ---------------------------------------------------------------------
hmri_config         = cfg_exbranch;
hmri_config.tag     = 'hmri_config';
hmri_config.name    = 'Configure toolbox';
hmri_config.val     = { hmri_setdef };
hmri_config.help    = {['Customised default parameters can be set here by selecting ' ...
    'a customised [hmri_local_defaults_*.m] file. Type [help hmri_local_defaults] for ' ...
    'more details.']};
hmri_config.prog    = @hmri_run_config;

end
%----------------------------------------------------------------------

% =========================================================================
% (VOUT &) RUN SUBFUNCTION(S)
% =========================================================================
function out = hmri_run_config(job)
%==========================================================================
% PURPOSE
% Load standard defaults and overwrite thme by customised values if
% provided via local defaults file. 
%==========================================================================

% reinitialise standard defaults
hmri_defaults;
% overwrite with customised defaults
if isfield(job.hmri_setdef,'customised')
    deffnam = job.hmri_setdef.customised;
    spm('Run',deffnam);
end
out = hmri_get_defaults;

end
%_______________________________________________________________________
