function create_mpm = tbx_scfg_hmri_create
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with the creation of the maps
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips

%--------------------------------------------------------------------------
% Imperfect spoiling correction - no correction
%--------------------------------------------------------------------------
iscnone         = cfg_entry;
iscnone.tag     = 'iscnone';
iscnone.name    = 'None';
iscnone.help    = {'No correction will be applied.'};
iscnone.strtype = 's';
iscnone.num     = [1 Inf];
iscnone.val     = {'none'};

%--------------------------------------------------------------------------
% Imperfect spoiling correction - select coefficients file
%--------------------------------------------------------------------------
iscfile         = cfg_files;
iscfile.tag     = 'iscfile';
iscfile.name    = 'Select ISC file';
iscfile.help    = {['Select an *.m file containing the correction coefficients ' ...
    'corresponding to the current acquisition protocol.']};
iscfile.filter = 'm';
iscfile.dir     = fullfile(fileparts(mfilename('fullpath')),'config');
iscfile.ufilter = '.*';
iscfile.num     = [1 1];

%--------------------------------------------------------------------------
% Imperfect spoiling correction - according to Phase Increment
%--------------------------------------------------------------------------
iscphase         = cfg_entry;
iscphase.tag     = 'iscphase';
iscphase.name    = 'Specify phase increment';
iscphase.help    = {['Specify the phase increment of the spoiled ' ...
     ' gradient echo sequence to correct imperfect spoiling. '...
     ' Possible values (in degree): 50 (Siemens), 117 (General Electrics), and 150 (Philips).']};
%scphase.values     = {50 117 150};
iscphase.val     = {50};

%--------------------------------------------------------------------------
% Imperfect spoiling correction - select coefficients
%--------------------------------------------------------------------------
isc        = cfg_choice;
isc.tag    = 'isc';
isc.name   = 'Imperfect spoiling correction (ISC)';
isc.help   = {['The RF- and gradient-spoiled gradient echo sequences used to ' ...
    'acquire the multiparametric mapping (MPM) data apply RF and gradient spoiling to ' ...
    'destroy unwanted transverse magnetisation. Imperfect spoiling can leave ' ...
    'residual bias in the apparent R1 map if no further correction is applied ' ...
    '[Preibisch and Deichmann 2009]. Correction coefficients are sequence-' ...
    'specific and can be determined by simulation [Preibisch and Deichmann 2009] to ' ...
    'account for deviation from the Ernst equation. '],[''], ...
    'The following options are available:', ...
    '- None [default]: no correction applied. ', ...
    ['- Select ISC file: you must choose a set of correction coefficients. ' ...
    'For the most standard MPM protocols using customised FLASH sequences on Siemens ' ...
    'scanners, the hMRI-toolbox provides spoiling correction coefficients. ' ...
    'For other protocols, the coefficients can be calculated seperately ' ...
    'and loaded here.'], ...
    '- Specify RF spoiling phase increment for correction according to Baudrexel et al. 2018',[''], ...
    ['WARNING: although the TR ' ...
    'and FA values are key parameters in simulating the imperfect spoiling, ' ...
    'correction coefficients cannot be chosen based on TR and FA only. ' ...
    'Amplitude and duration of gradient spoilers, phase increment for RF ' ...
    'spoiling and diffusion effects must be taken into account. ' ...
    'Therefore, two different sequences using identical TR anf FA are unlikely to ' ...
    'use identical correction coefficients.'],[''], ...
    ['ISC-MODULE: an additional hMRI module that efficiently calculates the ' ...
    'protocol-specific correction parameters required to account for ' ...
    'imperfect spoiling is planned. To be continued...']};
isc.values = {iscnone iscfile iscphase};
isc.val = {iscnone};


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
    'the Matlab Command Window, in a log file saved in the "Results/Supplementary" ' ...
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


% ---------------------------------------------------------------------
% Input FLASH images - T1-weighted 
% ---------------------------------------------------------------------
raws3           = cfg_files;
raws3.tag       = 'T1';
raws3.name      = 'T1 images';
raws3.help      = {'Input T1-weighted images.'};
raws3.filter    = 'image';
raws3.ufilter   = '.*';
raws3.num       = [0 Inf];
raws3.val       = {''};
% ---------------------------------------------------------------------
% Input FLASH images - PD-weighted
% ---------------------------------------------------------------------
raws2           = cfg_files;
raws2.tag       = 'PD';
raws2.name      = 'PD images';
raws2.help      = {'Input PD-weighted images.'};
raws2.filter    = 'image';
raws2.ufilter   = '.*';
raws2.num       = [0 Inf];
raws2.val       = {''};
% ---------------------------------------------------------------------
% Input FLASH images - MT-weighted
% ---------------------------------------------------------------------
raws1           = cfg_files;
raws1.tag       = 'MT';
raws1.name      = 'MT images';
raws1.help      = {'Input MT-weighted images.'};
raws1.filter    = 'image';
raws1.ufilter   = '.*';
raws1.num       = [0 Inf];
raws1.val       = {''};
% ---------------------------------------------------------------------
% All multiparameter input images
% ---------------------------------------------------------------------
raws            = cfg_branch;
raws.tag        = 'raw_mpm';
raws.name       = 'Multiparameter input images';
raws.help       = {'Input all the MT/PD/T1-weighted images.'};
raws.val        = {raws1 raws2 raws3};


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
% UNICORT B1 bias correction
% ---------------------------------------------------------------------
b1_input_UNICORT           = cfg_branch;
b1_input_UNICORT.tag       = 'UNICORT';
b1_input_UNICORT.name      = 'UNICORT';
b1_input_UNICORT.help      = {'UNICORT will be applied for B1 bias correction.'
    'No B1 input data required.'
    ['Customized processing parameters may be introduced by loading ' ...
    'customized defaults from a selected [hmri_b1_local_defaults_*.m] file.']};
b1_input_UNICORT.val       = {b1parameters};

% ---------------------------------------------------------------------
% No B1 bias correction
% ---------------------------------------------------------------------
b1_input_noB1           = cfg_entry;
b1_input_noB1.tag       = 'no_B1_correction';
b1_input_noB1.name      = 'no B1 correction';
b1_input_noB1.help      = {'No B1 bias correction will be applied.'
    ['NOTE: when no B1 map is available, UNICORT might be a better ' ...
    'solution than no B1 bias correction at all.']};
b1_input_noB1.strtype = 's';
b1_input_noB1.num     = [1 Inf];
b1_input_noB1.val     = {'noB1'};

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
% SDAM B1 protocol
% ---------------------------------------------------------------------
b1_input_SDAM           = cfg_branch;
b1_input_SDAM.tag       = 'SDAM';
b1_input_SDAM.name      = 'Saturated Double Angle Method';
b1_input_SDAM.help      = {'Saturated Double Angle Method (SDAM) protocol.', ...
    'As B1 input, please select a 2*alpha/alpha (e.g. 120°/60°) pair of images in that order.', ...
    ['Regarding processing parameters, you can either stick with metadata and standard ' ...
    'defaults parameters (recommended) or select your own [hmri_b1_local_defaults_*.m] customised defaults file ' ...
    '(fallback for situations where no metadata are available).']};
b1_input_SDAM.val       = {b1raw b1parameters};


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
    'toolbox when creating the multiparameter maps. See list below for a ' ...
    'brief description of each type. Note that all types may not be ' ...
    'available at your site.']
    [' - 3D EPI: B1map obtained from spin echo (SE) and stimulated echo ' ...
    '(STE) images recorded with a 3D EPI scheme [Lutti A et al., ' ...
    'PLoS One 2012;7(3):e32379].']
    [' - 3D AFI: 3D actual flip angle imaging (AFI) method based on [Yarnykh VL, ' ...
    'Magn Reson Med 2007;57:192-200].']
    [' - Saturated Double Angle Method (SDAM)']
    [' - tfl_b1_map: Siemens product sequence for B1 mapping based on turbo FLASH.']
    [' - rf_map: Siemens product sequence for B1 mapping based on SE/STE.']
    [' - no B1 correction: if selected no B1 bias correction will be applied.']
    [' - pre-processed B1: B1 map pre-calculated outside the hMRI toolbox, must ' ...
    'be expressed in percent units of the nominal flip angle value (percent bias).']
    [' - UNICORT: Use this option when B1 maps not available. ' ...
    'Bias field estimation and correction will be performed ' ...
    'using the approach described in [Weiskopf et al., NeuroImage 2011; 54:2116-2124]. ' ...
    'WARNING: the correction only applies to R1 maps.']
    }; %#ok<*NBRAK>
b1_type.values  = {b1_input_3DEPI b1_input_3DAFI b1_input_SDAM b1_input_tfl b1_input_rfmap b1_input_preproc b1_input_UNICORT b1_input_noB1};
b1_type.val     = {b1_input_3DEPI};

% ---------------------------------------------------------------------
% Input images for RF sensitivity - RF sensitivity maps for MTw images
% ---------------------------------------------------------------------
sraws3MT          = cfg_files;
sraws3MT.tag      = 'raw_sens_MT';
sraws3MT.name     = 'RF sensitivity maps for MTw images';
sraws3MT.help     = {'Select low resolution RF sensitivity maps acquired with the head and body coils respectively, in that order.'};
sraws3MT.filter   = 'image';
sraws3MT.ufilter  = '.*';
% sraws3MT.num      = [2 2];
sraws3MT.num       = [0 2];
sraws3MT.val       = {''};
% ---------------------------------------------------------------------
% Input images for RF sensitivity - RF sensitivity maps for PDw images
% ---------------------------------------------------------------------
sraws3PD          = cfg_files;
sraws3PD.tag      = 'raw_sens_PD';
sraws3PD.name     = 'RF sensitivity maps for PDw images';
sraws3PD.help     = {'Select low resolution RF sensitivity maps acquired with the head and body coils respectively, in that order.'};
sraws3PD.filter   = 'image';
sraws3PD.ufilter  = '.*';
% sraws3PD.num      = [2 2];
sraws3PD.num       = [0 2];
sraws3PD.val       = {''};
% ---------------------------------------------------------------------
% Input images for RF sensitivity - RF sensitivity maps for T1w images
% ---------------------------------------------------------------------
sraws3T1          = cfg_files;
sraws3T1.tag      = 'raw_sens_T1';
sraws3T1.name     = 'RF sensitivity maps for T1w images';
sraws3T1.help     = {'Select low resolution RF sensitivity maps acquired with the head and body coils respectively, in that order.'};
sraws3T1.filter   = 'image';
sraws3T1.ufilter  = '.*';
% sraws3T1.num      = [2 2];
sraws3T1.num       = [0 2];
sraws3T1.val       = {''};
% ---------------------------------------------------------------------
% xNULL No RF sensitivity bias correction applied at all
% ---------------------------------------------------------------------
xNULL         = cfg_entry;
xNULL.tag     = 'RF_none';
xNULL.name    = 'None';
xNULL.help    = {'No RF sensitivity bias correction will be applied.'};
xNULL.strtype = 's';
xNULL.num     = [1 Inf];
xNULL.val     = {'-'};
% ---------------------------------------------------------------------
% x0 No RF sensitivity
% ---------------------------------------------------------------------
x0         = cfg_entry;
x0.tag     = 'RF_us';
x0.name    = 'Unified Segmentation';
x0.help    = {['RF sensitivity bias correction based on the Unified Segmentation ' ...
    '(US) approach. The resulting Bias Field estimate is used to correct for ' ...
    'RF sensitivity bias (applies to the PD map calculation only). ' ...
    'No RF sensitivity map is required.']};
x0.strtype = 's';
x0.num     = [1 Inf];
x0.val     = {'-'};
% ---------------------------------------------------------------------
% x1 Single RF sensitivity maps acquired for all contrasts
% ---------------------------------------------------------------------
x1          = cfg_files;
x1.tag      = 'RF_once';
x1.name     = 'Single';
x1.help     = {'Single set of RF sensitivity maps acquired for all contrasts.', ...
    'Select low resolution RF sensitivity maps acquired with the head and body coils respectively, in that order.'};
x1.filter   = 'image';
x1.ufilter  = '.*';
x1.num      = [2 2];
% ---------------------------------------------------------------------
% x3 RF sensitivity acquired for each modality 
% ---------------------------------------------------------------------
x3         = cfg_branch;
x3.tag     = 'RF_per_contrast';
x3.name    = 'Per contrast';
x3.help    = {['One set of RF sensitivity maps is acquired for each contrast ' ...
    'i.e. for each of the PD-, T1- and MT-weighted multi-echo FLASH acquisitions.']};
x3.val     = {sraws3MT sraws3PD sraws3T1};
% ---------------------------------------------------------------------
% sensitivity Sensitivity choice
% ---------------------------------------------------------------------
sensitivity         = cfg_choice;
sensitivity.tag     = 'sensitivity';
sensitivity.name    = 'RF sensitivity bias correction';
sensitivity.help    = {'Specify the type of RF sensitivity bias correction to be applied. '
    'You can select either:'
    '- None: no correction will be applied,'
    '- Unified Segmentation: based on US, no RF sensitivity map required,'
    '- Single: based on a single set of RF sensitivity maps for all contrasts,'
    '- Per contrast: based on one set of RF sensitivity maps acquired for each contrast.'};
sensitivity.values  = {xNULL x0 x1 x3};
sensitivity.val     = {x0};

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
output.values  = {indir outdir };
output.val = {indir};

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj            = cfg_branch;
subj.tag        = 'subj';
subj.name       = 'Subject';
subj.help       = {'Specify a subject for maps calculation.'};
subj.val        = {output sensitivity b1_type raws isc popup};

% ---------------------------------------------------------------------
% data Data
% ---------------------------------------------------------------------
sdata           = cfg_repeat;
sdata.tag       = 'data';
sdata.name      = 'Few Subjects';
sdata.help      = {'Specify the number of subjects.'};
sdata.num       = [1 Inf];
sdata.val       = {subj };
sdata.values    = {subj };

% ---------------------------------------------------------------------
% create_mpr Create MPR maps (whether B0/B1 maps are available or not)
% ---------------------------------------------------------------------
create_mpm         = cfg_exbranch;
create_mpm.tag     = 'create_mpm';
create_mpm.name    = 'Create hMRI maps';
create_mpm.val     = { sdata };
create_mpm.help    = {'hMRI map creation based on multi-echo FLASH sequences including optional receive/transmit bias correction.'};
create_mpm.prog    = @hmri_run_create;
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
cdep(1,5*numel(job.subj)) = cfg_dep;
for i=1:numel(job.subj)
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('R1_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','R1','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('R2s_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','R2s','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('MT_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','MT','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('A_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','A','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('T1w_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','T1w','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('MTw_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','MTw','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('PDw_subj%d',i);
    cdep(k).src_output = substruct('.','subj','()',{i},'.','PDw','()',{':'});
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    k=k+1;
    
end
dep = cdep;
    
end
%_______________________________________________________________________
