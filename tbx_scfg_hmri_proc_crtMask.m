function proc_crtMask = tbx_scfg_hmri_proc_crtMask
% Configuration file for the mask creation routine
%_______________________________________________________________________
% Copyright (C) 2019 Cyclotron Research Center

% Written by Christophe Phillips

% ---------------------------------------------------------------------
% vols_smwc Smooth modulated warped tissue segement volumes
% ---------------------------------------------------------------------
vols_smwc         = cfg_files;
vols_smwc.tag     = 'vols_smwc';
vols_smwc.name    = 'smwc images';
vols_smwc.help    = {'Select the smooth modulated warped tissue segements (smwc*).', ... 
    'Pick only one type of smwc* images across all subjects!.', ...
    'Using ''spm_select'' recursive selection feature can help here.'};
vols_smwc.filter  = 'image';
vols_smwc.ufilter = '^smwc.*';
vols_smwc.num     = [1 Inf];

% ---------------------------------------------------------------------
% m_SMWCs Smooth modulate warped tissue segement (SMWC) maps
% ---------------------------------------------------------------------
m_SMWCs            = cfg_repeat;
m_SMWCs.tag        = 'm_SMWCs';
m_SMWCs.name       = 'Modulated warped tissue segements';
m_SMWCs.values     = {vols_smwc };
m_SMWCs.val        = {vols_smwc };
m_SMWCs.num = [1 Inf];
m_SMWCs.help       = {['Select the smooth modulated warped tissue segments ',...
    'of interest from all subjects.'], ...
    ['For the typical case of GM and WM, you would selectall the smwc1* images ', ...
    'in one set of ''smwc_images'' and the smwc2* ones in second set of ', ...
    '''smwc_images''!']};

% ---------------------------------------------------------------------
% threshTC Threshold for the mean smoothed modulated warped TCs
% ---------------------------------------------------------------------
threshTC         = cfg_entry;
threshTC.tag     = 'threshTC';
threshTC.name    = 'Mean smwc* threshold';
threshTC.val     = {.2};
threshTC.strtype = 'e';
threshTC.num     = [1 1];
threshTC.help    = {['Specify the threshold for the mean smoothed ' ...
    'modulated warped tissue class maps (smwc* images).']};


%--------------------------------------------------------------------------
% noOverlap Objective Function
%--------------------------------------------------------------------------
noOverlap         = cfg_menu;
noOverlap.tag     = 'noOverlap';
noOverlap.name    = 'Prevent mask overlap';
noOverlap.help    = {
    ['In order to prevent tissue mask from overlapping, i.e. having voxels ', ...
    'in both the GM and WM mask, each voxel can be assisgned to the ', ...
    'tissue class with the largest value.']
    }';
noOverlap.labels  = {
                    'Yes'
                    'No'
}';
noOverlap.values  = {
                    true
                    false
}';
noOverlap.val     = {true};

% ---------------------------------------------------------------------
% indir Input directory of 1st smwc image as output directory
% ---------------------------------------------------------------------
indir         = cfg_menu;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {['Output files will be written to the same folder as ',...
    'the 1st smwc image.']};
indir.labels  = {'Yes'};
indir.values  = {1};
indir.val     = {1};

% ---------------------------------------------------------------------
% outdir Output directory for all resulting images
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory, all together';
outdir.help    = {['Select a directory where all output files '... 
    'will be written to.']};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

% ---------------------------------------------------------------------
% output Output choice
% ---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {['Output directory can be the same as that of the 1st ',...
    'smwc image file or user selected.']};
output.values  = {indir outdir };
output.val     = {indir};

%--------------------------------------------------------------------------
% options Mask creation options
%--------------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Mask creation options';
options.val     = {threshTC noOverlap output};
options.help    = {''};

%% EXEC function
% ---------------------------------------------------------------------
% proc_crtMask Processing hMRI -> create mask image(s)
% ---------------------------------------------------------------------
proc_crtMask         = cfg_exbranch;
proc_crtMask.tag     = 'proc_crtMask';
proc_crtMask.name    = 'Proc. hMRI -> Mask creation';
proc_crtMask.val     = { m_SMWCs options };
proc_crtMask.check   = @check_proc_crtMask;
proc_crtMask.help    = { 
    ['Creating tissue specific mask images for further explicit masking ', ...
    'in the SPM analysis.']};
proc_crtMask.prog = @hmri_run_proc_crtMask;
proc_crtMask.vout = @vout_crtMask;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%__________________________________________________________________________
% Collect and prepare output
function dep = vout_crtMask(job)
% There should be one series of mask images and one series of mean images,
% 1 each per tissue class in input

n_TCs = numel(job.vols_smwc);     % #tissue classes

cdep = cfg_dep;
for ii=1:n_TCs
    cdep(end+1) = cfg_dep; %#ok<*AGROW>
    cdep(end).sname = sprintf('mask_TC #%d', ii);
    cdep(end).src_output = substruct('.', 'fn_maskTC', '()', {ii,':'});
    cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
end
for ii=1:n_TCs
    cdep(end+1) = cfg_dep;
    cdep(end).sname = sprintf('mean_TC #%d', ii);
    cdep(end).src_output = substruct('.', 'fn_meanTC', '()', {ii,':'});
    cdep(end).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
end
dep = cdep(2:end);

end

%__________________________________________________________________________
% Check number of files match
function chk = check_proc_crtMask(job)
% ensure they are the same for each list of files, one of each per subject.

n_TCs = numel(job.vols_smwc);
chk = '';

if n_TCs>1
    ni_TC = numel(job.vols_smwc{1});
    for ii=2:n_TCs
        if ni_TC ~= numel(job.vols_smwc{ii})
            chk = 'Incompatible number of tissue segments.';
        end
    end
end

end

