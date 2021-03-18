function create_mpm = tbx_scfg_hmri_B1_create
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with the creation of the maps
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips

%--------------------------------------------------------------------------
% To enable/disable pop-up messages for all warnings - recommended when
% piloting the data processing.
%--------------------------------------------------------------------------
popup        = cfg_menu;
popup.tag    = 'popup';
popup.name   = 'Pop-up warnings';
popup.help   = {['The user can review and keep track of all the information ' ...
    'collected, including warnings and other messages coming up during ' ...
    'the creation of the maps. By default, the information is logged in ' ...
    'the Matlab Command Window, in a log file saved in the "Results" ' ...
    'directory, and when more critical, displayed as a pop-up message.'], ...
    ['The latter must be disabled for processing series of datasets (since it ' ...
    'blocks the execution of the code) but it is strongly recommended to ' ...
    'leave it enabled when piloting the data processing (single subject) ' ...
    'to read through and acknowledge every message and make sure ' ...
    'everything is set up properly before running the processing on a ' ...
    'whole group.'], ...
    ['More information about the various messages and action to be taken ' ...
    '(or not) accordingly can be found on the hMRI-Toolbox WIKI (http://hmri.info). ' ...
    'In particular, see the "Debug tips & tricks" section.']};
popup.labels = {'Disable' 'Enable'};
popup.values = {false true};
popup.val = {true};

%--------------------------------------------------------------------------
% B1 acq/proc defaults file
%--------------------------------------------------------------------------
b1defaults         = cfg_files;
b1defaults.tag     = 'b1defaults';
b1defaults.name    = 'Customised B1 defaults file';
b1defaults.help    = {['Select the [hmri_b1_local_defaults_*.m] file containing ' ...
    'the parameters to process the B1 map data. By default, parameters will be ' ...
    'collected from metadata when available. Defaults parameters are provided as ' ...
    'fallback solution when metadata are not available and/or uncomplete.'], ...
    ['Please make sure that the parameters defined in the defaults file ' ...
    'are correct for your data. To create your own customised defaults file, ' ...
    'edit the distributed version and save it with a meaningful name such as ' ...
    '[hmri_b1_local_defaults_*myprotocol*.m].']};
b1defaults.filter  = 'm';
b1defaults.dir     = fullfile(fileparts(mfilename('fullpath')),'config','local');
b1defaults.ufilter = '^hmri_.*\.m$';
b1defaults.num     = [1 1];

% ---------------------------------------------------------------------
% Use metadata or standard defaults (no customization)
% ---------------------------------------------------------------------
b1metadata           = cfg_entry;
b1metadata.tag       = 'b1metadata';
b1metadata.name      = 'Use metadata or standard defaults';
b1metadata.help      = {''};
b1metadata.strtype = 's';
b1metadata.num     = [1 Inf];
b1metadata.val     = {'yes'};

%--------------------------------------------------------------------------
% B1 processing parameters
%--------------------------------------------------------------------------
b1parameters         = cfg_choice;
b1parameters.tag     = 'b1parameters';
b1parameters.name    = 'Processing parameters';
b1parameters.help    = {['You can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1parameters.values  = {b1metadata b1defaults};
b1parameters.val     = {b1metadata};

% ---------------------------------------------------------------------
% B1 input images 
% ---------------------------------------------------------------------
b1raw          = cfg_files;
b1raw.tag      = 'b1input';
b1raw.name     = 'B1 input';
b1raw.help     = {'Select B1 input images according to the type of B1 bias correction.'};
b1raw.filter   = 'image';
b1raw.ufilter  = '.*';
b1raw.num      = [2 30];
% b1raw.val      = {''};

% ---------------------------------------------------------------------
% B0 input images
% ---------------------------------------------------------------------
b0raw          = cfg_files;
b0raw.tag      = 'b0input';
b0raw.name     = 'B0 input';
b0raw.help     = {'Select B0 field map input images.' ...
    'Only required for distortion correction of EPI-based B1 maps.' ...
    'Select both magnitude images and the presubtracted phase image, in that order.'};
b0raw.filter   = 'image';
b0raw.ufilter  = '.*';
b0raw.num      = [3 3];
% b0raw.num      = [0 3];
% b0raw.val      = {''};

% ---------------------------------------------------------------------
% pre-calculated B1 map - including potential rescaling factor
% ---------------------------------------------------------------------
scafac         = cfg_entry;
scafac.tag     = 'scafac';
scafac.name    = 'Scaling factor';
scafac.help    = {'The values in the input B1 map will be multiplied by the provided factor.', ...
    ['If the input B1 map is already in percent units (p.u.) of the nominal flip angle ' ...
    'no need to apply any extra scaling factor (ScaFac = 1). If the input B1 map is a multiplication factor ' ...
    'of the nominal flip angle (i.e. value of 1 corresponds to the nominal flip angle), a ' ...
    'scaling factor ScaFac = 100 is required to produce a B1 map in p.u. of the ' ...
    'nominal flip angle.']};
scafac.strtype = 'r';
scafac.num     = [1 1];
scafac.val     = {1};

b1_input_preproc           = cfg_branch;
b1_input_preproc.tag       = 'pre_processed_B1';
b1_input_preproc.name      = 'pre-processed B1';
b1_input_preproc.help      = {'Input pre-calculated B1 bias map.'
    ['Please select one unprocessed magnitude image ' ...
    'from the B1map data set (for coregistration with the multiparameter maps) ' ...
    'and the preprocessed B1map, in that order.']
    ['The B1 map is expected to be in ' ...
    'percent units (p.u.) of the nominal flip angle. If this is not the case, ' ...
    'a scaling factor can be introduced (see Scaling factor description for more details).']};
b1_input_preproc.val       = {b1raw scafac};

% ---------------------------------------------------------------------
% RF_MAP B1 protocol
% ---------------------------------------------------------------------
b1_input_rfmap           = cfg_branch;
b1_input_rfmap.tag       = 'rf_map';
b1_input_rfmap.name      = 'rf_map';
b1_input_rfmap.help      = {'Input B1 images for rf_map B1 map protocol.' ...
    'As B1 input, please select the pair of anatomical and precalculated B1 map, in that order.'};
b1_input_rfmap.val       = {b1raw};

% ---------------------------------------------------------------------
% TFL_B1_MAP B1 protocol
% ---------------------------------------------------------------------
b1_input_tfl           = cfg_branch;
b1_input_tfl.tag       = 'tfl_b1_map';
b1_input_tfl.name      = 'tfl_b1_map';
b1_input_tfl.help      = {'Input B1 images for TFL B1 map protocol.' ...
    'As B1 input, please select the pair of anatomical and precalculated B1 map, in that order.'};
b1_input_tfl.val       = {b1raw};

% ---------------------------------------------------------------------
% i3D_AFI B1 protocol
% ---------------------------------------------------------------------
b1_input_3DAFI           = cfg_branch;
b1_input_3DAFI.tag       = 'i3D_AFI';
b1_input_3DAFI.name      = '3D AFI';
b1_input_3DAFI.help      = {'3D Actual Flip Angle Imaging (AFI) protocol.', ...
    'As B1 input, please select a TR2/TR1 pair of magnitude images.', ...
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_b1_local_defaults_*.m] customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1_input_3DAFI.val       = {b1raw b1parameters};

% ---------------------------------------------------------------------
% i3D_EPI B1 protocol
% ---------------------------------------------------------------------
b1_input_3DEPI           = cfg_branch;
b1_input_3DEPI.tag       = 'i3D_EPI';
b1_input_3DEPI.name      = '3D EPI';
b1_input_3DEPI.help      = {'Input B0/B1 data for 3D EPI protocol'
    'As B1 input, please select all pairs of SE/STE 3D EPI images.'
    ['For this EPI protocol, it is recommended to acquire B0 field map data ' ...
    'for distortion correction. If no B0 map available, the script will proceed ' ...
    'with distorted images.']
    ['Please enter the two magnitude images and the presubtracted phase image ' ...
    'from the B0 mapping acquisition, in that order.']
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_b1_local_defaults_*.m] customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1_input_3DEPI.val       = {b1raw b0raw b1parameters};

% ---------------------------------------------------------------------
% menu type_b1
% ---------------------------------------------------------------------
b1_type         = cfg_choice;
b1_type.tag     = 'b1_type';
b1_type.name    = 'B1 bias correction';
b1_type.help    = {'Choose the methods for B1 bias correction.'
    ['Various types of B1 mapping protocols can be handled by the hMRI ' ...
    'toolbox. See list below for a ' ...
    'brief description of each type. Note that all types may not be ' ...
    'available at your site.']
    [' - 3D EPI: B1map obtained from spin echo (SE) and stimulated echo ' ...
    '(STE) images recorded with a 3D EPI scheme [Lutti A et al., ' ...
    'PLoS One 2012;7(3):e32379].']
    [' - 3D AFI: 3D actual flip angle imaging (AFI) method based on [Yarnykh VL, ' ...
    'Magn Reson Med 2007;57:192-200].']
    [' - tfl_b1_map: Siemens product sequence for B1 mapping based on turbo FLASH.']
    [' - rf_map: Siemens product sequence for B1 mapping based on SE/STE.']
    [' - pre-processed B1: B1 map pre-calculated outside the hMRI toolbox, must ' ...
    'be expressed in percent units of the nominal flip angle value (percent bias).']
    }; %#ok<*NBRAK>
b1_type.values  = {b1_input_3DEPI b1_input_3DAFI b1_input_tfl b1_input_rfmap b1_input_preproc};
b1_type.val     = {b1_input_3DEPI};

% ---------------------------------------------------------------------
% indir Input directory as output directory
% ---------------------------------------------------------------------
indir         = cfg_entry;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {['Output files will be written to the same folder ' ...
    'as each corresponding input file.']};
indir.strtype = 's';
indir.num     = [1 Inf];
indir.val     = {'yes'};
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where output files will be written to.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% output Output choice
% ---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {['Output directory can be the same as the input ' ...
    'directory for each input file or user selected']};
output.values  = {indir outdir};
output.val = {indir};

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj            = cfg_branch;
subj.tag        = 'subj';
subj.name       = 'Subject';
subj.help       = {'Specify a subject for maps calculation.'};
subj.val        = {output b1_type popup};

% ---------------------------------------------------------------------
% data Data
% ---------------------------------------------------------------------
sdata           = cfg_repeat;
sdata.tag       = 'data';
sdata.name      = 'Few Subjects';
sdata.help      = {'Specify the number of subjects.'};
sdata.num       = [1 Inf];
sdata.val       = {subj};
sdata.values    = {subj};

% ---------------------------------------------------------------------
% create_mpr Create MPR maps (whether B0/B1 maps are available or not)
% ---------------------------------------------------------------------
create_mpm         = cfg_exbranch;
create_mpm.tag     = 'create_B1';
create_mpm.name    = 'Create B1 map';
create_mpm.val     = { sdata };
create_mpm.help    = {'Transmit bias field map creation for several common acquisition methods.'};
create_mpm.prog    = @hmri_run_create_B1;
create_mpm.vout    = @vout_create;

end
%----------------------------------------------------------------------

% ========================================================================
%% VOUT & OTHER SUBFUNCTIONS
% ========================================================================
% The RUN function:
% - out = hmri_run_create(job)
% is defined separately.
%_______________________________________________________________________

function dep = vout_create(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

k=1;
cdep(1,2*numel(job.subj)) = cfg_dep;
for i=1:numel(job.subj)
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('B1ref_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','B1ref','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('B1map_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','B1map','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
end
dep = cdep;
    
end
%_______________________________________________________________________