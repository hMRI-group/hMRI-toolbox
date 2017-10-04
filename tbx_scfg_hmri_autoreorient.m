function autoreor = tbx_scfg_hmri_autoreorient
% Configuration file for the "histological MRI" (hMRI) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% 
% PURPOSE: Reorientation of the images towards the MNI space is a standard
% step in neuroimage processing, and often a prerequisite for successful
% segmentation. We provide you with a simple tool for reorientation of all
% images prior to any further processing (including multiparameter map
% calculation). 
%
% METHODS: Reorientation is based on rigid-body coregistration of a
% suitable image (i.e. contrast must be well enough defined to allow for
% reliable coregistration) and application of the coregistration matrix to
% all images acquired during the same session. The code makes use of
% spm_affreg and templates available in SPM.
%
% AUTHORS: 
% Written by Christophe Phillips, 2011.
% Modified and extended by Evelyne Balteau for the hMRI Toolbox, 2017.
% Cyclotron Research Centre, University of Liege, Belgium

% ---------------------------------------------------------------------
% image Image to reorient
% ---------------------------------------------------------------------
ref         = cfg_files;
ref.tag     = 'reference';
ref.name    = 'Image';
ref.help    = {['Select the image that is best suited for rigid-body ' ...
    'coregistration to MNI space. Ideally, the GM/WM contrast in that image ' ...
    'should be well defined.']}';
ref.filter = 'image';
ref.ufilter = '.*';
ref.num     = [1 1];

% ---------------------------------------------------------------------
% imgtype Objective Function
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template';
template.help    = {['Auto-Reorient needs to know which modality of image ' ...
    'is selected and which template to use for the rigid-body coregistration.'] ,...
    ['You can select the template image from SPM12/canonical, SPM12/toolbox/OldNorm or any ' ...
    'suitable source of your choice.'], ...
    'By defaults, SPM12/canonical/avg152T1.nii is used.', ...
    'NOTE: The template must already be in the MNI orientation to make the reorientation effective!'};
template.filter  = 'image';
template.dir     = fullfile(spm('dir'),'canonical');
template.num     = [1 1];
template.def     = @(val)hmri_get_defaults('autoreorient_template', val{:});

%--------------------------------------------------------------------------
% other Other Images
%--------------------------------------------------------------------------
other         = cfg_files;
other.tag     = 'other';
other.name    = 'Other Images';
other.val     = {{''}};
other.help    = {['Select all the other images that need to remain in alignment ' ...
    'with the reference image that is reoriented. Typically, all the images acquired ' ...
    'during a single MRI session should be reoriented together. For convenience, select ' ...
    'all these images here, including the reference image.']};
other.filter  = 'image';
other.ufilter = '.*';
other.num     = [0 Inf];

%---------------------------------------------------------------------
% indir Input directory as output directory
%---------------------------------------------------------------------
indir         = cfg_entry;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {'Auto-reorientation is applied to the input files directly.'};
indir.strtype = 's';
indir.num     = [1 Inf];
indir.val     = {'yes'};
%---------------------------------------------------------------------
% outdir Output directory
%---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Input files are first copied to the selected directory before auto-reorientation.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
%---------------------------------------------------------------------
% output Output choice
% either "indir" (reorientation is applied to the input files directly) or
% user-defined "outdir" (images are copied to the "outdir" directory before
% reorientation). 
%---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {['Output directory can be the same as the input ' ...
    'directory or user selected.']};
output.values  = {indir outdir };
output.val = {indir};

%---------------------------------------------------------------------
% dep The way dependencies are made available
%---------------------------------------------------------------------
dep         = cfg_menu;
dep.tag     = 'dep';
dep.name    = 'Dependencies';
dep.help    = {['Dependencies can be made available as individual output ' ...
    'files or as a group of images, as is most convenient for the next ' ...
    'processing steps...']};
dep.labels = {'Individual', 'Grouped'}';
dep.values = {'individual', 'grouped'}';
dep.val = {'individual'};

%---------------------------------------------------------------------
% autoreor Auto-reorient one (or more) image(s) approximately in MNI space
%---------------------------------------------------------------------
autoreor         = cfg_exbranch;
autoreor.tag     = 'autoreor';
autoreor.name    = 'Auto-Reorient';
autoreor.val     = {ref template other output dep};
autoreor.help    = {[...
    'Function to automatically (but approximately) rigid-body reorient '...
    'a T1 image (or any other usual image modality) in the MNI space, '...
    'i.e. mainly set the AC location and correct for head rotation, in '...
    'order to further proceed with the segmentation/normalisation of the '...
    'image. This is useful since the Unified Segmentation process is '...
    'rather sensitive to the initial orientation of the image.'],'',...
    ['A set of other images can be reoriented along the 1st one. They '...
    'should be specified as "Other Images".'],'',...
    'The following outputs are available as dependencies for further processing:',...
    '- reoriented images (provided as individual images or as a group of images),',...
    ['- transformation matrix M (if needed to be applied to other images ' ...
    'at a later stage using SPM>Util>Reorient Images),'],...
    ['- inverted transformation matrix M (to go back to the initial '...
    'orientation using SPM>Util>Reorient Images).'],'',...
    ['NOTE: the job structure and the output structure are saved as '...
    'JSON data along with the reoriented images.']
    }';
autoreor.prog = @hmri_run_autoreorient;
autoreor.vout = @vout_autoreorient;
end

%------------------------------------------------------------------------
function dep = vout_autoreorient(job) %#ok<*INUSD>
% two options:
% - each reoriented image can be made available as a seperate dependency,
%   so it can be next "plugged" into the map creation module, or...
% - all reoriented images are grouped as a single dependency.

if strcmp(job.dep,'grouped')
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Auto-reoriented Image(s)';
    dep(1).src_output = substruct('.','files');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    cdep = 1;
else
    dep(numel(job.other)+1,1) = cfg_dep;
    % first thing first...
    fnam = spm_file(char(job.reference),'basename');
    dep(1)            = cfg_dep;
    dep(1).sname      = fnam;
    dep(1).src_output = substruct('.','files','{}',{1});%,'()',{':'});
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
    % then the others...
    for cdep=2:numel(job.other)+1
        fnam = spm_file(char(job.other{cdep-1}),'basename');
        dep(cdep)            = cfg_dep;
        dep(cdep).sname      = fnam;
        dep(cdep).src_output = substruct('.','files','{}',{cdep});%,'()',{':'});
        dep(cdep).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

% adding the transformation matrix M and its inverse invM to the
% dependencies
cdep = cdep+1;
dep(cdep)            = cfg_dep;
dep(cdep).sname      = 'transformation matrix M';
dep(cdep).src_output = substruct('.','M');

cdep = cdep+1;
dep(cdep)            = cfg_dep;
dep(cdep).sname      = 'inverted transformation matrix inv(M)';
dep(cdep).src_output = substruct('.','invM');

end
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function out = hmri_run_autoreorient(job)

ref = spm_file(char(job.reference),'number','');
template = spm_file(char(job.template),'number','');
other = spm_file(char(job.other),'number','');
output = job.output;
Nother = 0;
if ~isempty(other); Nother = size(other,1); end

% create AutoReorient directory and copy input files to it if required
if isfield(output,'outdir')
    outpath = fullfile(output.outdir{1},'AutoReorient'); % case outdir
    if ~exist(outpath,'dir'); mkdir(outpath); end
    copyref = fullfile(outpath, spm_file(ref,'filename'));
    copyfile(ref,copyref);
    try copyfile([spm_str_manip(ref,'r') '.json'],[spm_str_manip(copyref,'r') '.json']); end %#ok<*TRYNC>
    ref = copyref;
    copyother = cell(Nother,1);
    for cother=1:Nother
        copyother{cother} = fullfile(outpath, spm_file(other(cother,:),'filename'));
        copyfile(other(cother,:),copyother{cother});
        try copyfile([spm_str_manip(other(cother,:),'r') '.json'],[spm_str_manip(copyother{cother},'r') '.json']); end %#ok<*TRYNC>
    end
    other = char(copyother);
else 
    outpath = fileparts(ref);
end

out = hmri_autoreorient(ref, template, other);

if strcmp(job.dep,'individual')
    for cout=1:length(out.files)
        out.files{cout} = cellstr(out.files{cout});
    end
end

spm_jsonwrite(fullfile(outpath,'AutoReorient_job.json'),job,struct('indent','\t'));
spm_jsonwrite(fullfile(outpath,'AutoReorient_output.json'),out,struct('indent','\t'));

return;
end
%------------------------------------------------------------------------
