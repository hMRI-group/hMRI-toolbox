function proc_dart = tbx_scfg_hmri_proc_Dartel
% Configuration file for the "histological MRI" (hMRI) toolbox
% Applying DARTEL segmentation on previously segmented images
% There are 3 sub-modules:
% - Run Dartel (create Templates)
% - Run Dartel (use existing Templates)
% - Normalise to MNI
% The organisation of the _cfg_ and _run_ are very close to the orginal
% files from SPM (see Dartel toolbox itself) but are also tailored for the
% hMRI toolbox, i.e. the processing of a set of parametric maps from a
% series of subject.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips
% but largely inspired by the batch from the past VBQ toolbox.

% ---------------------------------------------------------------------
% Extract whole Dartel tbx configuration
% ---------------------------------------------------------------------
cfg_dartel = tbx_cfg_dartel;

% ---------------------------------------------------------------------
% indir Input directory as output directory
% ---------------------------------------------------------------------
indir         = cfg_menu;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {['Output files will be written to the same folder as ',...
    'each corresponding input file.']};
indir.labels  = {'Yes'};
indir.values  = {1};
indir.val     = {1};

% ---------------------------------------------------------------------
% outdir Output directory for all data
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory, all together';
outdir.help    = {['Select a directory where all output files from all '... 
    'subjects put together will be written to.']};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

% ---------------------------------------------------------------------
% outdir_ps Output directory for per-subject organisation
% ---------------------------------------------------------------------
outdir_ps         = cfg_files;
outdir_ps.tag     = 'outdir_ps';
outdir_ps.name    = 'Output directory, with per-subject sub-directory';
outdir_ps.help    = {['Select a directory where output files will be '...
    'written to, in each subject''s sub-directory.']};
outdir_ps.filter = 'dir';
outdir_ps.ufilter = '.*';
outdir_ps.num     = [1 1];

% ---------------------------------------------------------------------
% output Output choice
% ---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {['Output directory can be the same as the input ',...
    'directory for each input file or user selected (one for everything ',...
    'or preserving a per-subject organisation).']};
output.values  = {indir outdir outdir_ps };
output.val     = {indir};

% ---------------------------------------------------------------------
%% warp_dartel_cr Run DARTEL (create Templates) -> extract from whole config
% ---------------------------------------------------------------------
eval(['warp_dartel_cr = cfg_dartel', ...
    cfg_expr_values(cfg_dartel, 'warp'),';']);

% Add the output directory choice to the job structure, in 2nd position
% -> move the 'settings' to 3rd position
nFields = numel(warp_dartel_cr.val); %#ok<*NODEF>
if nFields>1
    for ii = nFields:-1:2
        warp_dartel_cr.val{ii+1} = warp_dartel_cr.val{ii};
    end
end
output_dartel_cr = output;

% % Adjust 'output' for Dartel_Create -> set the help of each option
% help_txt = {...
%     'Template files will be written to the same folder as 1st subject input file.', ...
%     'Select a directory where the Templates files will be saved.', ...
%     'Select a directory in which a sub-directory ''Dartel_Templates'' will contain the Templates files.'};
% name_txt = {...
%     'Same folder as 1st subject input file', ...
%     'Output directory, all together', ...
%     'Output directory, with ''DartelTemplates'' sub-directory'};
% for ii=1:3
%     output_dartel_cr.values{ii}.name = name_txt{ii};
%     output_dartel_cr.values{ii}.help = help_txt(ii);
% end
warp_dartel_cr.val{2} = output_dartel_cr;
warp_dartel_cr.prog  = @hmri_run_proc_dartel_template;

% ---------------------------------------------------------------------
%% warp_dartel_ex Run DARTEL (existing Templates) -> extract from whole config
% ---------------------------------------------------------------------
eval(['warp_dartel_ex = cfg_dartel', ...
    cfg_expr_values(cfg_dartel, 'warp1'),';']);

% Add the output directory choice to the job structure, in 2nd position
% -> move the 'settings' to 3rd position
nFields = numel(warp_dartel_ex.val); %#ok<*NODEF>
if nFields>1
    for ii = nFields:-1:2
        warp_dartel_ex.val{ii+1} = warp_dartel_ex.val{ii};
    end
end
output_dartel_cr = output;

warp_dartel_ex.val{2} = output_dartel_cr;
warp_dartel_ex.prog  = @hmri_run_proc_dartel_warp;

% ---------------------------------------------------------------------
% vols_field Deformation fields
% ---------------------------------------------------------------------
vols_field         = cfg_files;
vols_field.tag     = 'vols_field';
vols_field.name    = 'Flow fields';
vols_field.help    = {'Flow fields.'};
vols_field.filter  = 'image';
vols_field.ufilter = '^u_rc.*';
vols_field.num     = [1 Inf];

% ---------------------------------------------------------------------
% vols_pm Parametric volumes
% ---------------------------------------------------------------------
vols_pm         = cfg_files;
vols_pm.tag     = 'vols_pm';
vols_pm.name    = 'Maps';
vols_pm.help    = {['Select whole brain parameter maps (e.g. MT, R2*, ',...
    'FA etc), one type from all subjects (keep the order consistent!).']};
vols_pm.filter  = 'image';
vols_pm.ufilter = '.*';
vols_pm.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_pams Parameter maps sets
% ---------------------------------------------------------------------
m_pams        = cfg_repeat;
m_pams.tag    = 'm_pams';
m_pams.name   = 'Parameter maps';
m_pams.values = {vols_pm };
m_pams.val    = {vols_pm };
m_pams.num    = [1 Inf];
m_pams.help   = {['Select whole brain parameter maps (e.g. MT, ',...
    'R2*, FA etc).']};

% ---------------------------------------------------------------------
% vols_tc Tissue segments
% ---------------------------------------------------------------------
vols_tc         = cfg_files;
vols_tc.tag     = 'vols_tc';
vols_tc.name    = 'c* images';
vols_tc.help    = {'Select the tissue segment images (c*)'};
vols_tc.filter  = 'image';
vols_tc.ufilter = '^c[\d].*'; % filenames starting with 'c' and a number
vols_tc.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_TCs Tissue segment image sets
% ---------------------------------------------------------------------
m_TCs            = cfg_repeat;
m_TCs.tag        = 'maps';
m_TCs.name       = 'Segmented tissue class';
m_TCs.values     = {vols_tc };
m_TCs.val        = {vols_tc };
m_TCs.num = [0 Inf];
m_TCs.help       = {['Select the tissue segement images',...
    'of interest from all subjects. This should typically be the c1* ',...
    'and c2* images, for GM and WM, in 2 separate sets of images.']};

% ---------------------------------------------------------------------
% multsdata Data
% ---------------------------------------------------------------------
multsdata = cfg_branch;
multsdata.tag = 'multsdata';
multsdata.name = 'Data';
multsdata.val = {m_TCs m_pams vols_field};
% multsdata.val = {multsdata_gm multsdata_wm multsdata_f multsdata_u};

% ---------------------------------------------------------------------
%% nrm Normalize to MNI -> extract from whole config
% ---------------------------------------------------------------------
eval(['nrm = cfg_dartel', ...
    cfg_expr_values(cfg_dartel, 'mni_norm'),';']);

% ---------------------------------------------------------------------
% Reorganize SPM-batch entries for the hMRI-batch
% ---------------------------------------------------------------------
eval(['nrm' cfg_expr(nrm, 'data') '= multsdata;']);
% Drop out 2 fields:
% - the 'preserve' field, as this is fixed in the run function
% - the 'fwhm' field, as no smoothin applied here
eval(['nrm' regexprep(cfg_expr(nrm, 'preserve'), '{([0-9]+)}$', '($1)') '=[];']);
eval(['nrm' regexprep(cfg_expr(nrm, 'fwhm'), '{([0-9]+)}$', '($1)') '=[];']);

% Add the output directory chooice to the job structure, in 3rd position
nFields = numel(nrm.val);
if nFields>2
    for ii = nFields:-1:3
        nrm.val{ii+1} = nrm.val{ii};
    end
end
nrm.val{3} = output;

% Function calls
nrm.prog  = @hmri_run_proc_dartel_norm;
nrm.vout  = @vout_norm_fun;
nrm.check = [];

% ---------------------------------------------------------------------
% dartel DARTEL Tools
% ---------------------------------------------------------------------
proc_dart         = cfg_choice;
proc_dart.tag     = 'proc_dart';
proc_dart.name    = 'Proc. hMRI -> Dartel';
proc_dart.values  = {warp_dartel_cr warp_dartel_ex nrm };
%dartel.num     = [0 Inf];

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

% =======================================================================
%% VOUT & CHECK FUNCTIONS
% =======================================================================

%_______________________________________________________________________

function dep = vout_norm_fun(job) %#ok<*INUSD>
cdep = cfg_dep;
% deal with mwc's
for ii=1:numel(job.multsdata.vols_tc)
    cdep(end+1) = cfg_dep;
    cdep(end).sname      = sprintf('mwc#%d',ii);
    cdep(end).src_output = substruct('.','vols_mwc','{}',{ii});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
% deal with wc's
for ii=1:numel(job.multsdata.vols_tc)
    cdep(end+1) = cfg_dep;
    cdep(end).sname      = sprintf('mwc#%d',ii);
    cdep(end).src_output = substruct('.','vols_wc','{}',{ii});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%deal with wPM's
for ii=1:numel(job.multsdata.vols_pm)
    cdep(end+1) = cfg_dep;
    cdep(end).sname      = sprintf('warped Param Map #%d',ii);
    cdep(end).src_output = substruct('.','vols_wpm','{}',{ii});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
dep = cdep(2:end);

end
%_______________________________________________________________________

function dep = vout_dartel_warp(job) %#ok<*DEFNU>
for i=1:numel(job.images{1})
    fdep(i)            = cfg_dep;
    fdep(i).sname      = sprintf('Flow Field_subj%d',i);
    fdep(i).src_output = substruct('.','files','()',{i});
    fdep(i).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
dep = fdep;
end
%_______________________________________________________________________

function dep = vout_dartel_template(job)

if isa(job.settings.template,'cfg_dep') || ~ ...
        isempty(deblank(job.settings.template))
    for it=0:numel(job.settings.param),
        tdep(it+1)            = cfg_dep;
        tdep(it+1).sname      = sprintf('Template (Iteration %d)', it);
        tdep(it+1).src_output = substruct('.','template','()',{it+1});
        tdep(it+1).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
else
    tdep = cfg_dep;
    tdep = tdep(false);
end

for i=1:numel(job.images{1})
    fdep(i)            = cfg_dep;
    fdep(i).sname      = sprintf('Flow Field_subj%d',i);
    fdep(i).src_output = substruct('.','files','()',{i});
    fdep(i).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
dep = [tdep fdep];
end
%_______________________________________________________________________

function chk = check_dartel_template(job)
n1 = numel(job.images);
n2 = numel(job.images{1});
chk = '';
for i=1:n1,
    if numel(job.images{i}) ~= n2,
        chk = 'Incompatible number of images';
        break;
    end;
end;
end

% ========================================================================
%% SUBFUNCTIONS to handle matlabbatch structure and fields
% ========================================================================

function expr = cfg_expr_values(c, varargin) %#ok<INUSL>
% Extracting the 'values' field with index matching the input argument and
% 'tag' of the matlabbatch structure.

expr = 'c';
for i=1:size(varargin,2)
    %         if strcmp(class(varargin{i}), 'double')
    if isa(varargin{i}, 'double')
        expr = [expr '.values{' num2str(varargin{i}) '}'];  %#ok<*AGROW>
    else
        v = eval([expr ';']);
        for j=1:size(v.values,2)
            if strcmp(v.values{j}.tag, varargin{i})
                break
            end
        end
        expr = [expr '.values{' num2str(j) '}']; 
    end
end
expr = expr(2:end);
end
%_______________________________________________________________________

function expr = cfg_expr(c, varargin) %#ok<INUSL>
expr = 'c';
for i=1:size(varargin,2)
    %         if strcmp(class(varargin{i}), 'double')
    if isa(varargin{i}, 'double')
        expr = [expr '.val{' num2str(varargin{i}) '}']; 
    else
        v = eval([expr ';']);
        for j=1:size(v.val,2)
            if strcmp(v.val{j}.tag, varargin{i})
                break
            end
        end
        expr = [expr '.val{' num2str(j) '}']; 
    end
end
expr = expr(2:end);
end
%_______________________________________________________________________

function c = unlimit(c)
try
    if isa(c, 'cfg_files')
        c.num = [0 Inf];
    end
catch e %#ok<*NASGU>
end
try
    for i=1:numel(c.val)
        c.val{i} = unlimit(c.val{i});
    end
catch e
end
end
%_______________________________________________________________________

function c = cfg_set_val(c, varargin)
expr = ['c' cfg_expr(c, varargin{1:end-1})];
eval([expr '.val={varargin{end}};']);
end
%_______________________________________________________________________

