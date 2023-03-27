function create_mpm = tbx_scfg_hmri_create
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

% ---------------------------------------------------------------------
% menu b1_type
% ---------------------------------------------------------------------
[b1_type,b1parameters] = tbx_scfg_hmri_B1_menu;

% ---------------------------------------------------------------------
% New B1 mapping methods should be added to tbx_scfg_hmri_B1_menu.m
% Only B1 mapping methods which cannot be run independently of the rest of 
% the hMRI toolbox should be added here!
% ---------------------------------------------------------------------

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
% Add extra B1 mapping methods which cannot be run independently of the
% MPM map creation to the menu
% ---------------------------------------------------------------------
b1_type.values  = [b1_type.values {b1_input_UNICORT, b1_input_noB1}];
b1_type.help=[b1_type.help; {...
    [' - UNICORT: Use this option when B1 maps not available. ' ...
    'Bias field estimation and correction will be performed ' ...
    'using the approach described in [Weiskopf et al., NeuroImage 2011; 54:2116-2124]. ' ...
    'WARNING: the correction only applies to R1 maps.']
    [' - no B1 correction: This option is *not* recommended when computing R1 maps. ' ...
    'Consider using UNICORT instead.']
    }];

% ---------------------------------------------------------------------
% Input images for RF sensitivity - RF sensitivity maps for MTw images
% ---------------------------------------------------------------------
sraws3MT          = cfg_files;
sraws3MT.tag      = 'raw_sens_MT';
sraws3MT.name     = 'RF sensitivity maps for MTw images';
sraws3MT.help     = {'Select low resolution calibration images.  ' ...
    'The first image should have been acquired immediately prior to the MTw acquisition.' ...
    'The second image is the reference, e.g. acquired with the body coil (see Papp et al. MRM 2016) immediately prior to the MTw acquisition, ' ...
    'or use the array coil image for a particular contrast (the same in each entry) as the reference (see Balbastre et al. MRM 2022).'};
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
sraws3PD.help     = {'Select low resolution calibration images.  ' ...
    'The first image should have been acquired immediately prior to the PDw acquisition.' ...
    'The second image is the reference, e.g. acquired with the body coil (see Papp et al. MRM 2016) immediately prior to the PDw acquisition, ' ...
    'or use the array coil image for a particular contrast (the same in each entry) as the reference (see Balbastre et al. MRM 2022).'};
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
sraws3T1.help     = {'Select low resolution calibration images.  ' ...
    'The first image should have been acquired immediately prior to the T1w acquisition.' ...
    'The second image is the reference, e.g. acquired with the body coil (see Papp et al. MRM 2016) immediately prior to the T1w acquisition, ' ...
    'or use the array coil image for a particular contrast (the same in each entry) as the reference (see Balbastre et al. MRM 2022).'};
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
subj.val        = {output sensitivity b1_type raws popup};

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