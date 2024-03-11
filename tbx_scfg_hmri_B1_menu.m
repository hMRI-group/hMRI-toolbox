function [b1_type,b1parameters] = tbx_scfg_hmri_B1_menu
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with menu options for the creation of B1 maps
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips


%--------------------------------------------------------------------------
% B1 acq/proc defaults file
%--------------------------------------------------------------------------
b1defaults         = cfg_files;
b1defaults.tag     = 'b1defaults';
b1defaults.name    = 'Customised B1 defaults file';
b1defaults.help    = {['Select the [hmri_b1_*_defaults_*.m] file containing ' ...
    'the parameters to process the B1 map data. By default, parameters will be ' ...
    'collected from metadata when available. Defaults parameters are provided as ' ...
    'fallback solution when metadata are not available and/or uncomplete.'], ...
    ['Please make sure that the parameters defined in the defaults file ' ...
    'are correct for your data. To create your own customised defaults file, ' ...
    'edit the distributed version and save it with a meaningful name such as ' ...
    '[hmri_b1_local_defaults_*myprotocol*.m].']};
b1defaults.filter  = 'm';
b1defaults.dir     = fullfile(fileparts(mfilename('fullpath')),'config','local');
b1defaults.ufilter = [...
    '(?i)',   ... ignore case
    '^hmri_', ... filename starts with 'hmri_'
    '.*b1.*', ... filename includes 'b1' somewhere in the middle
    '\.m$'    ... filename ends with '.m'
    ];
b1defaults.num     = [1 1];

% ---------------------------------------------------------------------
% Use metadata or standard defaults (no customization)
% ---------------------------------------------------------------------
b1metadata         = cfg_entry;
b1metadata.tag     = 'b1metadata';
b1metadata.name    = 'Use metadata or standard defaults';
b1metadata.help    = {''};
b1metadata.strtype = 's';
b1metadata.num     = [1 Inf];
b1metadata.val     = {'yes'};

%--------------------------------------------------------------------------
% B1 processing parameters
%--------------------------------------------------------------------------
b1parameters        = cfg_choice;
b1parameters.tag    = 'b1parameters';
b1parameters.name   = 'Processing parameters';
b1parameters.help   = {['You can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1parameters.values = {b1metadata b1defaults};
b1parameters.val    = {b1metadata};

% ---------------------------------------------------------------------
% B1 input images 
% ---------------------------------------------------------------------
% b1raw is the default covering most used cases. 
% If specific methods require changes to these defaults, please create 
% a new field with a unique name and use that.
b1raw         = cfg_files;
b1raw.tag     = 'b1input';
b1raw.name    = 'B1 input';
b1raw.help    = {'Select B1 input images according to the type of B1 bias correction.'};
b1raw.filter  = 'image';
b1raw.ufilter = '.*';
b1raw.num     = [2 2];

% 3DEPI can have more than 2 input images
b1raw_3DEPI     = b1raw;
b1raw_3DEPI.num = [2 30];

% ---------------------------------------------------------------------
% B0 input images
% ---------------------------------------------------------------------
b0raw         = cfg_files;
b0raw.tag     = 'b0input';
b0raw.name    = 'B0 input';
b0raw.help    = {'Select B0 field map input images.' ...
    'Only required for distortion correction of EPI-based B1 maps.' ...
    'Select both magnitude images and the presubtracted phase image, in that order.'};
b0raw.filter  = 'image';
b0raw.ufilter = '.*';
b0raw.num     = [3 3];

% ---------------------------------------------------------------------
% pre-calculated B1 map - including potential rescaling factor
% ---------------------------------------------------------------------
scafac         = cfg_entry;
scafac.tag     = 'scafac';
scafac.name    = 'Scaling factor';
scafac.help    = {'The values in the input B1 map will be multiplied by the provided factor.'
    ['1. If the input B1 map is already in percent units (p.u.) of the nominal flip angle ' ...
    'then there is no need to apply any extra scaling factor (set ScaFac = 1).'] 
    ['2. If the input B1 map is a multiplication factor ' ...
    'for the nominal flip angle (i.e. a value of 1 corresponds to the nominal flip angle), ' ...
    'then set ScaFac = 100 to produce a B1 map in p.u.']
    ['3. If the input map is a flip angle map with nominal flip angle alpha, then ' ...
    'set ScaFac = 100/alpha to produce a B1 map in p.u. of the nominal flip angle.']};
scafac.strtype = 'r';
scafac.num     = [1 1];
scafac.val     = {1};

b1_input_preproc      = cfg_branch;
b1_input_preproc.tag  = 'pre_processed_B1';
b1_input_preproc.name = 'pre-processed B1';
b1_input_preproc.help = {'Input pre-calculated B1 bias map.'
    ['Please select one unprocessed magnitude image ' ...
    'from the B1map data set (for coregistration with the multiparameter maps) ' ...
    'and the preprocessed B1map, in that order.']
    ['The B1 map is expected to be in ' ...
    'percent units (p.u.) of the nominal flip angle. If this is not the case, ' ...
    'a scaling factor can be introduced (see Scaling factor description for more details).']};
b1_input_preproc.val  = {b1raw scafac b1parameters};


% ---------------------------------------------------------------------
% RF_MAP B1 protocol
% ---------------------------------------------------------------------
b1_input_rfmap      = cfg_branch;
b1_input_rfmap.tag  = 'rf_map';
b1_input_rfmap.name = 'rf_map';
b1_input_rfmap.help = {'Input B1 images for rf_map B1 map protocol.' ...
    'As B1 input, please select the pair of anatomical and precalculated B1 map, in that order.'};
b1_input_rfmap.val  = {b1raw b1parameters};


% ---------------------------------------------------------------------
% TFL_B1_MAP B1 protocol
% ---------------------------------------------------------------------
b1_input_tfl      = cfg_branch;
b1_input_tfl.tag  = 'tfl_b1_map';
b1_input_tfl.name = 'tfl_b1_map';
b1_input_tfl.help = {'Input B1 images for TFL B1 map protocol.' ...
    'As B1 input, please select the pair of anatomical and precalculated B1 map, in that order.'};
b1_input_tfl.val  = {b1raw b1parameters};


% ---------------------------------------------------------------------
% SDAM B1 protocol
% ---------------------------------------------------------------------
b1_input_DAM      = cfg_branch;
b1_input_DAM.tag  = 'DAM';
b1_input_DAM.name = 'Double Angle Method';
b1_input_DAM.help = {'Double Angle Method (DAM) protocol.' ...
    'As B1 input, please select a 2*alpha/alpha (e.g. 120°/60°) pair of images in that order.' ...
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_b1_local_defaults_*.m] customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1_input_DAM.val  = {b1raw b1parameters};


% ---------------------------------------------------------------------
% i3D_AFI B1 protocol
% ---------------------------------------------------------------------
b1_input_3DAFI      = cfg_branch;
b1_input_3DAFI.tag  = 'i3D_AFI';
b1_input_3DAFI.name = '3D AFI';
b1_input_3DAFI.help = {'3D Actual Flip Angle Imaging (AFI) protocol.' ...
    'As B1 input, please select the pair of TR1 and TR2 magnitude images.' ...
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_b1_local_defaults_*.m] customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1_input_3DAFI.val  = {b1raw b1parameters};


% ---------------------------------------------------------------------
% i3D_EPI B1 protocol
% ---------------------------------------------------------------------
b1_input_3DEPI      = cfg_branch;
b1_input_3DEPI.tag  = 'i3D_EPI';
b1_input_3DEPI.name = '3D EPI';
b1_input_3DEPI.help = {'Input B0/B1 data for 3D EPI protocol'
    'As B1 input, please select all pairs of SE/STE 3D EPI images.'
    ['The images are expected to be in acquisition order. This corresponds to SE_beta1, STE_beta1, SE_beta2, STE_beta2, ' ...
    '... SE_betaN, STE_betaN, where beta1 to betaN are the different flip angles in the same order as in "b1acq.beta". ']
    ['For qMRI-BIDS compatible data, the BIDS metadata "FlipAngle" and "EchoTime" ' ...
    'can be used to sort the data into SE and STE if the appropriate flags ' ...
    '("b1validation.checkTEs" and "b1validation.useBidsFlipAngleField") are "true" ' ...
    'in the b1 defaults file.']
    ['For this EPI protocol, it is recommended to acquire B0 field map data ' ...
    'for distortion correction. If no B0 map available, the script will proceed ' ...
    'with distorted images.']
    ['Please enter the two magnitude images and the presubtracted phase image ' ...
    'from the B0 mapping acquisition, in that order.']
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_b1_local_defaults_*.m] customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1_input_3DEPI.val  = {b1raw_3DEPI b0raw b1parameters};


% ---------------------------------------------------------------------
% menu b1_type
% ---------------------------------------------------------------------
b1_type        = cfg_choice;
b1_type.tag    = 'b1_type';
b1_type.name   = 'B1 bias correction';
b1_type.help   = {'Choose the methods for B1 bias correction.'
    ['Various types of B1 mapping protocols can be handled by the hMRI ' ...
    'toolbox. See list below for a ' ...
    'brief description of each type. Note that all types may not be ' ...
    'available at your site.']
    [' - 3D EPI: B1map obtained from spin echo (SE) and stimulated echo ' ...
    '(STE) images recorded with a 3D EPI scheme [Lutti A et al., ' ...
    'PLoS One 2012;7(3):e32379].']
    [' - 3D AFI: 3D actual flip angle imaging (AFI) method based on [Yarnykh VL, ' ...
    'Magn Reson Med 2007;57:192-200].']
    ' - Double Angle Method (DAM)'
    ' - tfl_b1_map: Siemens product sequence for B1 mapping based on turbo FLASH.'
    ' - rf_map: Siemens product sequence for B1 mapping based on SE/STE.'
    [' - pre-processed B1: B1 map pre-calculated outside the hMRI toolbox, must ' ...
    'be proportional to the nominal flip angle value.']};
b1_type.values = {b1_input_3DEPI b1_input_3DAFI b1_input_DAM b1_input_tfl b1_input_rfmap b1_input_preproc};
b1_type.val    = {b1_input_3DEPI};

end
