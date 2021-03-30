function create_B1 = tbx_scfg_hmri_B1_create
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with the creation of B1 maps
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

% ---------------------------------------------------------------------
% menu type_b1
% ---------------------------------------------------------------------
b1_type = tbx_scfg_hmri_B1_menu;

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
% create_B1 Create B1 maps
% ---------------------------------------------------------------------
create_B1         = cfg_exbranch;
create_B1.tag     = 'create_B1';
create_B1.name    = 'Create B1 map';
create_B1.val     = { sdata };
create_B1.help    = {'Transmit bias field map creation for several common acquisition methods.'};
create_B1.prog    = @hmri_run_create_B1;
create_B1.vout    = @vout_create;

end
%----------------------------------------------------------------------

% ========================================================================
%% VOUT & OTHER SUBFUNCTIONS
% ========================================================================
% The RUN function:
% - out = hmri_run_create_B1(job)
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