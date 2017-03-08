function proc_us = tbx_scfg_hmri_proc_US
% Configuration file for the segmentation part of the processing modules of
% the "histological MRI" (hMRI) toolbox.
% -> Apply "unifies segementation" (US) on series of images.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by Christophe Phillips
% but largely inspired by the batch from the past VBQ toolbox.

% NOTE:
% It could be advantageous to define the TPM in a definition file and use
% it when ever we need it. Right now, this is hard-coded in the cfg file!

% -------------------------------------------------------------------------
% vols Volumes
% ---------------------------------------------------------------------
vols            = cfg_files;
vols.tag        = 's_vols';
vols.name       = 'Reference (T1w or MT) images';
vols.help       = {'Select reference images, i.e. T1w or MT, for "unified segmentation".'};
vols.filter     = 'image';
vols.ufilter    = '.*';
vols.num        = [1 Inf];

% ---------------------------------------------------------------------
% vols_pm Parametric maps
% ---------------------------------------------------------------------
vols_pm         = cfg_files;
vols_pm.tag     = 'mp_vols';
vols_pm.name    = 'Volumes';
vols_pm.help    = {['Select whole brain parameter maps (e.g. MT, R2*, ',...
    'FA etc) for processing.']};
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
indir.labels = {'Yes'};
indir.values = {1};

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
output.help    = {['Output directory can be the same as the input ',...
    'directory for each input file or user selected']};
output.values  = {indir outdir };

% ---------------------------------------------------------------------
% Get stuff from SPM-US config file & adapting it
% ---------------------------------------------------------------------
% get the prproc8 config object
preproc8 = spm_cfg_preproc8;
% set the data->channel->vols to here defined 'vols'
eval(['preproc8',cfg_expr(preproc8, 'data', 'channel', 'vols'),' = vols;']);
% % set the output for first 3 tissue classes to write warped, mod & unmod.
% for i=1:3
%     eval(['preproc8',cfg_expr(preproc8, 'tissues', i, 'warped'),'.val{1}=[1 1];']);
% end

% ---------------------------------------------------------------------
% struct Structurals
% ---------------------------------------------------------------------
% extract the data channel from preproc8 for the structural reference def.
rstruct = eval(['preproc8',cfg_expr(preproc8, 'data', 'channel')]);
rstruct.tag = 'rstruct';
rstruct.name = 'Ref. structurals';
% set the bias cutoff to 'no bias' correction
rstruct = cfg_set_val(rstruct, 'biasfwhm', Inf);

% ---------------------------------------------------------------------
% many_pams Parameter maps
% ---------------------------------------------------------------------
% used for 'many subjects', i.e. list the data per map type across subjects
many_pams            = cfg_repeat;
many_pams.tag        = 'maps';
many_pams.name       = 'Parameter maps';
many_pams.values        = {vols_pm };
many_pams.val        = {vols_pm };
many_pams.num = [1 Inf];
many_pams.help       = {['Select whole brain parameter maps (e.g. MT, ',...
    'R2*, FA etc) for processing.']};

% ---------------------------------------------------------------------
% many_sdatas Many subjects data
% ---------------------------------------------------------------------
many_sdatas = cfg_branch;
many_sdatas.tag = 'many_sdatas';
many_sdatas.name = 'Many Subjects';
many_sdatas.val = {output many_pams unlimit(rstruct)};
many_sdatas.help = {'Specify images for many subjects'};

% ---------------------------------------------------------------------
% preproc8 Segment MT/T1w data
% ---------------------------------------------------------------------
proc_us = preproc8;
proc_us.name = 'Maps preprocessing - Segmentation';
proc_us.tag  = 'proc_us';
% Combine data defintion (local) with the tissue specs & other parameters 
% from preproc8 (these 2 are the last elements in preproc8.val)
proc_us.val = [{many_sdatas} preproc8.val(2:end)];

% set the output for the 6 tissue classes to
% - 1st 3 -> write warped, mod+unmod, and native, native+dartelImp. 
% - last 3 -> nothing
% plus update tpm with the hMRI specific
w_native = [[1 1];[1 1];[1 1];[0 0];[0 0];[0 0]];
w_warped = [[1 1];[1 1];[1 1];[0 0];[0 0];[0 0]];
fn_tpm = fullfile(spm('dir'),'toolbox','hMRI','tpm','unwTPM_sl2.nii');
for ii=1:size(w_native,1)
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'native', w_native(ii,:));
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'native', w_warped(ii,:));
    proc_us = cfg_set_val(proc_us, 'tissues', ii, 'tpm', ...
        {spm_file(fn_tpm,'number',ii)});
    
end

% set the output to write out the forward deformation field 
proc_us = cfg_set_val(proc_us, 'warp', 'write', [0 1]);

proc_us.prog = @hmri_run_local_preproc;
proc_us.vout = @vout_preproc;

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% =======================================================================
% VOUT function
% =======================================================================
% TO ADD:
% Need for a check function to ensure the same number of files in each
% series of maps + reference structural.

function dep = vout_preproc(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

cdep = cfg_dep;

% if isfield(job, 'many_few_sdatas')
%     if isfield(job.many_few_sdatas, 'subjc')
%         job.subjc = job.many_few_sdatas.subjc;
%     else
        for i=1:numel(job.tissue)
            if job.tissue(i).native(1)
                cdep(end+1) = cfg_dep;
                cdep(end).sname = sprintf('c%d Images', i);
                cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'c', '()', {':'});
                cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
            end
            if job.tissue(i).native(2)
                cdep(end+1) = cfg_dep;
                cdep(end).sname = sprintf('rc%d Images', i);
                cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'rc', '()', {':'});
                cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
            end
        end
        
%         disp(job.many_few_sdatas);
%         for i=1:numel(job.many_few_sdatas.many_sdatas.mp_vols)
        for i=1:numel(job.many_sdatas.mp_vols)
            cdep(end+1) = cfg_dep;
            cdep(end).sname = sprintf('%d Parameter Volumes', i);
            cdep(end).src_output = substruct('.', 'maps', '()', {i}, '.', 'mp_vols', '()', {':'});
            cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
        end
        
        dep = cdep(2:end);
end
%         return;
%     end
% end
% for nm=1:numel(job.subjc)
%     for i=1:numel(job.tissue),
%         if job.tissue(i).native(1),
%             cdep(end+1)          = cfg_dep; %#ok<*AGROW>
%             cdep(end).sname      = sprintf('c%d_subj%d Images',i,nm);
%             cdep(end).src_output = substruct('.','subjc','()',{nm},'.','tiss','()',{i},'.','c','()',{':'});
%             cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%         end
%         if job.tissue(i).native(2),
%             cdep(end+1)          = cfg_dep;
%             cdep(end).sname      = sprintf('rc%d_subj%d Images',i,nm);
%             cdep(end).src_output = substruct('.','subjc','()',{nm},'.','tiss','()',{i},'.','rc','()',{':'});
%             cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%         end
%     end
%     for i=1:numel(job.subjc(nm).maps.mp_vols)
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('%d_subj%d Parameter Volumes',i,nm);
%         cdep(end).src_output = substruct('.','subjc','()',{nm},'.','maps','.','mp_vols','()',{i});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
%     
% end
% 
% dep = cdep(2:end);
% 
% end


%% =======================================================================
% SUBFUNCTIONS to handle matlabbatch structure and fields
% ========================================================================

% function expr = cfg_expr_values(c, varargin) %#ok<INUSL>
% % Extracting the 'values' field with index matching the input argument and
% % 'tag' of the matlabbatch structure.
% 
% expr = 'c';
% for i=1:size(varargin,2)
%     %         if strcmp(class(varargin{i}), 'double')
%     if isa(varargin{i}, 'double')
%         expr = [expr '.values{' num2str(varargin{i}) '}']; 
%     else
%         v = eval([expr ';']);
%         for j=1:size(v.values,2)
%             if strcmp(v.values{j}.tag, varargin{i})
%                 break
%             end
%         end
%         expr = [expr '.values{' num2str(j) '}']; 
%     end
% end
% expr = expr(2:end);
% end
%_______________________________________________________________________

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

