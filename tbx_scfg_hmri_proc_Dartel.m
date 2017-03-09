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

% TO DO
% Define the dependencies and output from the "normalize to MNI module"!

% ---------------------------------------------------------------------
% warp Run DARTEL (create Templates)
% ---------------------------------------------------------------------
warp_dartel_cr = tbx_cfg_dartel;
eval(['warp_dartel_cr = warp_dartel_cr', ...
    cfg_expr_values(warp_dartel_cr, 'warp'),';']);

% ---------------------------------------------------------------------
% warp1 Run DARTEL (existing Templates)
% ---------------------------------------------------------------------
warp_dartel_ex = tbx_cfg_dartel;
eval(['warp_dartel_ex = warp_dartel_ex', ...
    cfg_expr_values(warp_dartel_ex, 'warp1'),';']);

% ---------------------------------------------------------------------
% vols_field Deformation fields
% ---------------------------------------------------------------------
vols_field         = cfg_files;
vols_field.tag     = 'vols_field';
vols_field.name    = 'Flow fields';
vols_field.help    = {'Flow fields.'};
vols_field.filter  = 'image';
vols_field.ufilter = '.*';
vols_field.num     = [1 Inf];

% ---------------------------------------------------------------------
% vols_pm Parametric volumes
% ---------------------------------------------------------------------
vols_pm         = cfg_files;
vols_pm.tag     = 'vols_mp';
vols_pm.name    = 'Volumes';
vols_pm.help    = {['Select whole brain parameter maps (e.g. MT, R2*, ',...
    'FA etc).']};
vols_pm.filter  = 'image';
vols_pm.ufilter = '.*';
vols_pm.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_pams Parameter maps, used for 'many subjects'
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
% vols_tc Parametric volumes
% ---------------------------------------------------------------------
vols_tc         = cfg_files;
vols_tc.tag     = 'vols_tc';
vols_tc.name    = 'rc* images';
vols_tc.help    = {'Select the Dartel imported tissue classes (rc*)'};
vols_tc.filter  = 'image';
vols_tc.ufilter = '^rc.*';
vols_tc.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_TCs Tissue class (TC) maps, used for 'many subjects'
% ---------------------------------------------------------------------
m_TCs            = cfg_repeat;
m_TCs.tag        = 'maps';
m_TCs.name       = 'Dartel imported tissue class';
m_TCs.values     = {vols_tc };
m_TCs.val        = {vols_tc };
m_TCs.num = [1 Inf];
m_TCs.help       = {['Select the modulated warped tissue classes (TC) ',...
    'of interest from all subjects. This ould typically be the mwc1* ',...
    'and mwc2* images for GM and WM.']};

% ---------------------------------------------------------------------
% multsdata Data
% ---------------------------------------------------------------------
multsdata = cfg_branch;
multsdata.tag = 'multsdata';
multsdata.name = 'Data';
multsdata.val = {m_TCs m_pams vols_field};
% multsdata.val = {multsdata_gm multsdata_wm multsdata_f multsdata_u};

nrm_tmp = tbx_cfg_dartel;
eval(['nrm = nrm_tmp', ...
    cfg_expr_values(nrm_tmp, 'mni_norm'),';']);

eval(['nrm' cfg_expr(nrm, 'data') '= multsdata;']); %#ok<NODEF>
% Drop out 2 fields:
% - the 'preserve' field, as this is fixed in the run function
% - the 'fwhm' field, as no smoothin applied here
eval(['nrm' regexprep(cfg_expr(nrm, 'preserve'), '{([0-9]+)}$', '($1)') '=[];']);
eval(['nrm' regexprep(cfg_expr(nrm, 'fwhm'), '{([0-9]+)}$', '($1)') '=[];']);

nrm.prog  = @hmri_run_local_dartel_norm_fun;
nrm.vout  = @vout_norm_fun;
nrm.check = [];

% ---------------------------------------------------------------------
% dartel DARTEL Tools
% ---------------------------------------------------------------------
proc_dart         = cfg_choice;
proc_dart.tag     = 'proc_dart';
proc_dart.name    = 'Proc. hMRI -> Dartel';
proc_dart.help    = {
    ['This toolbox is based around the ``A Fast Diffeomorphic ',...
    'Registration Algorithm'''' paper/* \cite{ashburner07} */. The idea ',...
    'is to register images by computing a ``flow field'''', which can ',...
    'then be ``exponentiated'''' to generate both forward and backward ',...
    'deformations. Currently, the software only works with images that ',...
    'have isotropic voxels, identical dimensions and which are in ',...
    'approximate alignment with each other. One of the reasons for this ',...
    'is that the approach assumes circulant boundary conditions, ',...
    'which makes modelling global rotations impossible. Another reason ',...
    'why the images should be approximately aligned is because there ',...
    'are interactions among the transformations that are minimised by ',...
    'beginning with images that are already almost in register. This ',...
    'problem could be alleviated by a time varying flow field, but ',...
    'this is currently computationally impractical.']
    ['Because of these limitations, images should first be imported. This ',...
    'involves taking the ``*_seg_sn.mat'''' files produced by the ',...
    'segmentation code of SPM12, and writing out rigidly transformed ',...
    'versions of the tissue class images, such that they are in as close ',...
    'alignment as possible with the tissue probability maps. Rigidly ',...
    'transformed original images can also be generated, with the option ',...
    'to have skull-stripped versions.']
    ['The next step is the registration itself.  This can involve ',...
    'matching single images together, or it can involve the simultaneous ',...
    'registration of e.g. GM with GM, WM with WM and 1-(GM+WM) with ',...
    '1-(GM+WM) (when needed, the 1-(GM+WM) class is generated ',...
    'implicitly, so there is no need to include this class yourself). ',...
    'This procedure begins by creating a mean of all the images, which ',...
    'is used as an initial template. Deformations from this template to ',...
    'each of the individual images are computed, and the template is ',...
    'then re-generated by applying the inverses of the deformations to ',...
    'the images and averaging. This procedure is repeated a number of times.',...
    'Finally, warped versions of the images (or other images that are ',...
    'in alignment with them) can be generated. ']
    ''
    ['This toolbox is not yet seamlessly integrated into the SPM package. ',...
    'Eventually, the plan is to use many of the ideas here as the ',...
    'default strategy for spatial normalisation. The toolbox may change ',...
    'with future updates.  There will also be a number of other (as ',...
    'yet unspecified) extensions, which may include a variable velocity ',...
    'version (related to LDDMM). Note that the Fast Diffeomorphism paper ',...
    'only describes a sum of squares objective function. The multinomial ',...
    'objective function is an extension, based on a more appropriate ',...
    'model for aligning binary data to a template.']
    }';
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
dep = cfg_dep;
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------------------------------------------------------------
% % configuration for STEP 1: Spatial preprocessing, i.e. segmentation
% % -------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % vols Volumes
% % ---------------------------------------------------------------------
% vols            = cfg_files;
% vols.tag        = 's_vols';
% vols.name       = 'T1w or MT images';
% vols.help       = {'Select T1w or MT images for "unified segmentation".'};
% vols.filter     = 'image';
% vols.ufilter    = '.*';
% vols.num        = [1 Inf];
% % ---------------------------------------------------------------------
% % Gaussian FWHM
% % ---------------------------------------------------------------------
% fwhm         = cfg_entry;
% fwhm.tag     = 'fwhm';
% fwhm.name    = 'Gaussian FWHM';
% fwhm.val     = {[6 6 6]};
% fwhm.strtype = 'e';
% fwhm.num     = [1 3];
% fwhm.help    = {'Specify the full-width at half maximum (FWHM) of the ',...
%     'Gaussian blurring kernel in mm. Three values should be entered',...
%     'denoting the FWHM in the x, y and z directions. Note that you can ',...
%     'also specify [0 0 0]',...
%     'but any ``modulated'' data will show aliasing (see eg Wikipedia), ',...
%     'which occurs because of the way the warped images are generated.'};
% % ---------------------------------------------------------------------
% % vols Volumes
% % ---------------------------------------------------------------------
% vols_pm         = cfg_files;
% vols_pm.tag     = 'mp_vols';
% vols_pm.name    = 'Volumes';
% vols_pm.help    = {'Select whole brain parameter maps (e.g. MT, R2*, ',...
%     'FA etc) for processing.'};
% vols_pm.filter  = 'image';
% vols_pm.ufilter = '.*';
% vols_pm.num     = [1 Inf];
% % ---------------------------------------------------------------------
% % pams Data = parametric maps
% % ---------------------------------------------------------------------
% pams            = cfg_branch;
% pams.tag        = 'maps';
% pams.name       = 'Parameter maps';
% pams.val        = {vols_pm };
% pams.help       = {'Select whole brain parameter maps (e.g. MT, R2*, ',...
%     'FA etc) for processing.'};
% % ---------------------------------------------------------------------
% % indir Input directory as output directory
% % ---------------------------------------------------------------------
% indir         = cfg_menu;
% indir.tag     = 'indir';
% indir.name    = 'Input directory';
% indir.help    = {'Output files will be written to the same folder as ',...
%     'each corresponding input file.'};
% indir.labels = {'Yes'};
% indir.values = {1};
% % ---------------------------------------------------------------------
% % outdir Output directory
% % ---------------------------------------------------------------------
% outdir         = cfg_files;
% outdir.tag     = 'outdir';
% outdir.name    = 'Output directory';
% outdir.help    = {'Select a directory where output files will be written to.'};
% outdir.filter = 'dir';
% outdir.ufilter = '.*';
% outdir.num     = [1 1];
% % ---------------------------------------------------------------------
% % output Output choice
% % ---------------------------------------------------------------------
% output         = cfg_choice;
% output.tag     = 'output';
% output.name    = 'Output choice';
% output.help    = {'Output directory can be the same as the input ',...
%     'directory for each input file or user selected'};
% output.values  = {indir outdir };
% % ---------------------------------------------------------------------
% % struct Structurals
% % ---------------------------------------------------------------------
% preproc8 = spm_cfg_preproc8;
% eval(['preproc8',cfg_expr(preproc8, 'data', 'channel', 'vols'),' = vols;']);
% for i=1:3
%     eval(['preproc8',cfg_expr(preproc8, 'tissues', i, 'warped'),'.val{1}=[1 1];']);
% end
% struct = eval(['preproc8',cfg_expr(preproc8, 'data', 'channel')]);
% struct.tag = 'struct';
% struct.name = 'Structurals';
% % ---------------------------------------------------------------------
% % subjc Subject, used for 'few subjects'
% % ---------------------------------------------------------------------
% subjc            = cfg_branch;
% subjc.tag        = 'subjc';
% subjc.name       = 'Subject';
% subjc.val        = {output pams struct };
% subjc.help       = {'Specify a subject for maps calculation.'};
% % ---------------------------------------------------------------------
% % sdatas Data for few subjects
% % ---------------------------------------------------------------------
% sdatas           = cfg_repeat;
% sdatas.tag       = 'data';
% sdatas.name      = 'Few subjects';
% sdatas.val       = {subjc };
% sdatas.help      = {'Specify the number of subjects. Note that all ',...
%     'raw images have to be entered in the order MT/PD/T1/B1/B0.'};
% sdatas.values    = {subjc };
% sdatas.num       = [1 Inf];
% % ---------------------------------------------------------------------
% % many_pams Parameter maps, used for 'many subjects'
% % ---------------------------------------------------------------------
% many_pams            = cfg_repeat;
% many_pams.tag        = 'maps';
% many_pams.name       = 'Parameter maps';
% many_pams.values        = {vols_pm };
% many_pams.val        = {vols_pm };
% many_pams.num = [1 Inf];
% many_pams.help       = {'Select whole brain parameter maps (e.g. MT, ',...
%     'R2*, FA etc) for processing.'};
% % ---------------------------------------------------------------------
% % many_sdatas Many subjects data
% % ---------------------------------------------------------------------
% many_sdatas = cfg_branch;
% many_sdatas.tag = 'many_sdatas';
% many_sdatas.name = 'Many Subjects';
% many_sdatas.val = {output many_pams unlimit(struct)};
% many_sdatas.help = {'Specify images for many subjects'};
% % ---------------------------------------------------------------------
% % Choice many/few
% % ---------------------------------------------------------------------
% many_few_sdatas = cfg_choice;
% many_few_sdatas.tag = 'many_few_sdatas';
% many_few_sdatas.name = 'Data Specification Method';
% many_few_sdatas.values = {sdatas, many_sdatas};
% many_few_sdatas.val = {sdatas};
% many_few_sdatas.help = {'Specify your data either as many or few subjects.'};
% % ---------------------------------------------------------------------
% % preproc8 Segment MT/T1w data
% % ---------------------------------------------------------------------
% preproc8.name = 'Maps preprocessing - Segmentation';
% preproc8.val = [{many_few_sdatas} preproc8.val(2:end) {fwhm}];
% 
% preproc8 = cfg_set_val(preproc8, 'tissues', 1, 'native', [1 1]);
% preproc8 = cfg_set_val(preproc8, 'tissues', 2, 'native', [1 1]);
% preproc8 = cfg_set_val(preproc8, 'tissues', 3, 'native', [1 1]);
% preproc8 = cfg_set_val(preproc8, 'tissues', 4, 'native', [0 0]);
% preproc8 = cfg_set_val(preproc8, 'tissues', 5, 'native', [0 0]);
% preproc8 = cfg_set_val(preproc8, 'tissues', 6, 'native', [0 0]);
% preproc8 = cfg_set_val(preproc8, 'warp', 'write', [0 1]);
% preproc8.prog = @hmri_run_local_preproc;
% preproc8.vout = @vout_preproc;



% % ---------------------------------------------------------------------
% % proc Preprocess maps
% % ---------------------------------------------------------------------
% proc_dart         = cfg_choice;
% proc_dart.tag     = 'proc_dart';
% proc_dart.name    = 'Proc. hMRI -> Dartel';
% proc_dart.help    = {
%     ['Parameter maps are registered to standard space, scaled and ready ',...
%     'for histological MRI (hMRI) analysis.']
%     }';
% proc_dart.values  = {proc_dart };
% 
% % ---------------------------------------------------------------------
% % multsdata_gm GM Images
% % ---------------------------------------------------------------------
% multsdata_gm         = cfg_files;
% multsdata_gm.tag     = 'multsdata_gm';
% multsdata_gm.name    = 'GM Volumes';
% multsdata_gm.help    = {'Select GM volumes.'};
% multsdata_gm.filter  = 'image';
% multsdata_gm.ufilter = '.*';
% multsdata_gm.num     = [1 Inf];
% % ---------------------------------------------------------------------
% % multsdata_wm WM Images
% % ---------------------------------------------------------------------
% multsdata_wm         = cfg_files;
% multsdata_wm.tag     = 'multsdata_wm';
% multsdata_wm.name    = 'WM Volumes';
% multsdata_wm.help    = {'Select WM volumes.'};
% multsdata_wm.filter  = 'image';
% multsdata_wm.ufilter = '.*';
% multsdata_wm.num     = [1 Inf];
% % ---------------------------------------------------------------------
% % multsdata_f1 Multi-parameter maps
% % ---------------------------------------------------------------------
% multsdata_f1         = cfg_files;
% multsdata_f1.tag     = 'multsdata_f1';
% multsdata_f1.name    = 'Map';
% multsdata_f1.help    = {'Select multi-parameter maps.'};
% multsdata_f1.filter  = 'image';
% multsdata_f1.ufilter = '.*';
% multsdata_f1.num     = [1 Inf];
% % ---------------------------------------------------------------------
% % multsdata_f1 Multi-parameter maps
% % ---------------------------------------------------------------------
% multsdata_f         = cfg_repeat;
% multsdata_f.tag     = 'multsdata_f';
% multsdata_f.name    = 'Multi-parameter maps';
% multsdata_f.val     = { multsdata_f1 };
% multsdata_f.help    = {'Select multi-parameter maps.'};
% multsdata_f.values = { multsdata_f1 };
% multsdata_f.num     = [1 Inf];


% function dep = vout_preproc(job)
% % This depends on job contents, which may not be present when virtual
% % outputs are calculated.
% 
% cdep = cfg_dep;
% 
% if isfield(job, 'many_few_sdatas')
%     if isfield(job.many_few_sdatas, 'subjc')
%         job.subjc = job.many_few_sdatas.subjc;
%     else
%         for i=1:numel(job.tissue)
%             if job.tissue(i).native(1)
%                 cdep(end+1) = cfg_dep;
%                 cdep(end).sname = sprintf('c%d Images', i);
%                 cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'c', '()', {':'});
%                 cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
%             end
%             if job.tissue(i).native(2)
%                 cdep(end+1) = cfg_dep;
%                 cdep(end).sname = sprintf('rc%d Images', i);
%                 cdep(end).src_output = substruct('.', 'tiss', '()', {i}, '.', 'rc', '()', {':'});
%                 cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
%             end
%         end
%         
%         disp(job.many_few_sdatas);
%         for i=1:numel(job.many_few_sdatas.many_sdatas.mp_vols)
%             cdep(end+1) = cfg_dep;
%             cdep(end).sname = sprintf('%d Parameter Volumes', i);
%             cdep(end).src_output = substruct('.', 'maps', '()', {i}, '.', 'mp_vols', '()', {':'});
%             cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
%         end
%         
%         dep = cdep(2:end);
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
