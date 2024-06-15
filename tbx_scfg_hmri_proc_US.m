function proc_us = tbx_scfg_hmri_proc_US
% Configuration file for the segmentation part of the processing modules of
% the "histological MRI" (hMRI) toolbox.
% -> Apply "unified segmentation" (US) on series of images.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips
% but largely inspired by the batch from the past VBQ toolbox.

% NOTE:
% It could be advantageous to define the TPM in a definition file and use
% it when ever we need it. Right now, this is hard-coded in the cfg file!

% ---------------------------------------------------------------------
% vols_pm Parametric maps
% ---------------------------------------------------------------------
vols_pm         = cfg_files;
vols_pm.tag     = 'vols_pm';
vols_pm.name    = 'Maps';
vols_pm.help    = {['Select whole brain maps of one type (e.g. MT, R2*, ',...
    'R1, etc.) from all subjects for processing.']};
vols_pm.filter  = 'image';
vols_pm.ufilter = '.*';
vols_pm.num     = [1 Inf];

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
    'or preserve a per-subject organisation).']};
output.values  = {indir outdir outdir_ps };
output.val     = {indir};

% ---------------------------------------------------------------------
% Get stuff from SPM-US config file & adapting it
% ---------------------------------------------------------------------
% get the prproc8 config object
preproc8 = spm_cfg_preproc8;
% implement Guillaume's fix 
% for compatibility with spm/spm (spm-development-version)
if isa(preproc8.val,'function_handle')
    preproc8.val = feval(preproc8.val);
end
% % set the bias cutoff and regularisation to 'no bias' 'no reg' correction
% preproc8 = cfg_set_val(preproc8, 'data', 'channel', 'biasfwhm', Inf); 
% preproc8 = cfg_set_val(preproc8, 'data', 'channel', 'biasreg', 0); 
% NOTE that this only change the values for the 1st channel!
% -> use dirty trick here under to both change 1st and default values...

% ---------------------------------------------------------------------
% rstruct Structurals
% ---------------------------------------------------------------------
% extract the data channel from preproc8 for the structural reference def.
rstruct = eval(['preproc8',cfg_expr(preproc8, 'data')]);
rstruct.tag = 'rstruct';
rstruct.name = 'Structurals for segmentation';

% set the bias cutoff to 'no bias' correction
rstruct = cfg_set_val(rstruct, 'channel', 'biasfwhm', Inf); 
% set the bias regularisation to 'no regularisation' correction
rstruct = cfg_set_val(rstruct, 'channel', 'biasreg', 0); 

% Dirty trick to set the default values of biasfhwm to Inf for new channels
kk_exp = cfg_expr(rstruct, 'channel', 'biasfwhm');
ll = strfind(kk_exp,'.val{1}'); % pick bits with '.val{1}'
% Replace 1nd occurence by '.values{1}' -> defaults
kk_exp = [kk_exp(1:ll(1)-1),'.values{1}',kk_exp((ll(1)+7):end)];
eval(['rstruct',kk_exp '.val = {Inf};']);
% Same goes for 'biasreg' sot to 0
kk_exp = cfg_expr(rstruct, 'channel', 'biasreg');
ll = strfind(kk_exp,'.val{1}'); % pick bits with '.val{1}'
% Replace 1nd occurence by '.values{1}' -> defaults
kk_exp = [kk_exp(1:ll(1)-1),'.values{1}',kk_exp((ll(1)+7):end)];
eval(['rstruct',kk_exp '.val = {0};']);

% ---------------------------------------------------------------------
% many_pams Parameter maps
% ---------------------------------------------------------------------
% used for 'many subjects', i.e. list the data per map type across subjects
many_pams            = cfg_repeat;
many_pams.tag        = 'maps';
many_pams.name       = 'Parametric maps';
many_pams.values     = {vols_pm};
many_pams.val        = {}; % Empty to begin with
many_pams.num = [0 Inf];
many_pams.help       = {['Select whole brain parameter maps (e.g. MT, ',...
    'R2*, FA, etc.) from all subjects for processing, one type per entry.']};

% ---------------------------------------------------------------------
% vox Voxel sizes
% ---------------------------------------------------------------------
vox          = cfg_entry;
vox.tag      = 'vox';
vox.name     = 'Voxel sizes';
vox.num      = [1 3];
vox.strtype  = 'e';
vox.val      = {[1 1 1]};
vox.help     = {[...
'Specify the voxel sizes of the deformation field and tissue classes ',...
'to be produced. Non-finite values will default to the voxel sizes of ',...
'the template image that was originally used to estimate the deformation.']};

%--------------------------------------------------------------------------
% bb Bounding box
%--------------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'};
bb.strtype = 'r';
bb.num     = [2 3];
bb.def     = @(val)spm_get_defaults('normalise.write.bb', val{:});

% ---------------------------------------------------------------------
% many_sdatas Many subjects data
% ---------------------------------------------------------------------
many_sdatas = cfg_branch;
many_sdatas.tag = 'many_sdatas';
many_sdatas.name = 'Data & options';
many_sdatas.val = {output unlimit(rstruct) many_pams vox bb};
many_sdatas.help = {'Specify images for many subjects at once.' , ...
    ['Processing will work on 1 subject at the time, using his ',...
    'structural image(s) to estimate the segmentation and warping parameters. ', ...
    'Then warps are applied *only* on his parametric maps, if provided.']};

% ---------------------------------------------------------------------
% preproc8 Segment MT/T1w data
% ---------------------------------------------------------------------
proc_us = preproc8;
proc_us.name = 'Proc. hMRI -> Segmentation';
proc_us.tag  = 'proc_us';
% Combine data defintion (local) with the tissue specs & other parameters
% from preproc8 (these 2 are the last elements in preproc8.val)
proc_us.val = [{many_sdatas} preproc8.val(2:end)];

% get the output for the tissue classes
w_native = hmri_get_defaults('proc.w_native');
w_warped = hmri_get_defaults('proc.w_warped');
% get the number of Gaussians per tissue class
nGauss = hmri_get_defaults('proc.nGauss');
% use the hMRI specific TPMs.
fn_tpm = hmri_get_defaults('proc.TPM');
% Fill in each tissue class parameters
for ii=1:size(w_native,1)
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'native', w_native(ii,:));
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'warped', w_warped(ii,:));
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'tpm', ...
        {spm_file(fn_tpm,'number',ii)});
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'ngaus', nGauss(ii));
end
% set the output to write out the forward deformation field
proc_us = cfg_set_val(proc_us, 'warp', 'write', [0 1]);

proc_us.prog = @hmri_run_proc_US;
proc_us.vout = @vout_preproc;
proc_us.check = @check_USdata;

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% =======================================================================
% VOUT function
% =======================================================================

function dep = vout_preproc(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

cdep = cfg_dep;

% Collect tissue class images (4 of them)
for i=1:numel(job.tissue)
    if job.tissue(i).native(1)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('c%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'c', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}}); % cfg_findspec({{'filter','nifti'}});
    end
    if job.tissue(i).native(2)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('rc%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'rc', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).warped(1)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('wc%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'wc', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tissue(i).warped(2)
        cdep(end+1) = cfg_dep;
        cdep(end).sname = sprintf('mwc%d Images', i);
        cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'mwc', '()', {':'});
        cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%         cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

% Collect warped parametric maps
for i=1:numel(job.many_sdatas.vols_pm)
    cdep(end+1) = cfg_dep;
    cdep(end).sname = sprintf('Warped p. maps #%d', i);
    cdep(end).src_output = substruct('.', 'maps', '()', {i}, '.', 'wvols_pm', '()', {':'});
    cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
%     cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
end

% Collect the deformation fields
cdep(end+1) = cfg_dep;
cdep(end).sname = 'Def. fields';
cdep(end).src_output = substruct('.', 'def', '.', 'fn', '()', {':'});
cdep(end).tgt_spec = cfg_findspec({{'filter','nifti'}});
% cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});

dep = cdep(2:end);
end

%% =======================================================================
% CHECKING the data
% ========================================================================
function t = check_USdata(job)
% Checking that the data are consistent.
t   = {};

nSubj = numel(job.many_sdatas.channel(1).vols);
nChan = numel(job.many_sdatas.channel); % number of channels
nPara = numel(job.many_sdatas.vols_pm); % number of maps type
% 1/ Check the number of structurals in each channel
if nChan>1
    for ii=2:nChan
        if numel(job.many_sdatas.channel(ii).vols)~=0 && ...
                numel(job.many_sdatas.channel(ii).vols)~=nSubj
        t{1} = 'Structural channels have different number of images/subjects!';
            warndlg(t,'Structural channel numbers');
            return
        end    
    end
end
% 2/ Check this number matches the number of parametric maps
if nPara>1
    for ii=1:nPara
        if numel(job.many_sdatas.vols_pm{ii})~=0 && ...
                numel(job.many_sdatas.vols_pm{ii})~=nSubj
        t{1} = 'Number of maps not matching number of structural images/subjects!';
            warndlg(t,'Maps numbers');
            return
        end    
    end
end

end

%% =======================================================================
% SUBFUNCTIONS to handle matlabbatch structure and fields
% ========================================================================

function expr = cfg_expr(c, varargin) %#ok<INUSL>
expr = 'c';
for i=1:size(varargin,2)
    %         if strcmp(class(varargin{i}), 'double')
    if isa(varargin{i}, 'double')
        expr = [expr '.val{' num2str(varargin{i}) '}'];  %#ok<*AGROW>
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

function c = cfg_set_val(c, varargin)
expr = ['c' cfg_expr(c, varargin{1:end-1})];
eval([expr '.val={varargin{end}};']);
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

