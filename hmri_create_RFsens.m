function jobsubj = hmri_create_RFsens(jobsubj)
% RF sensitivity calculation as part of the hmri toolbox
% Based on a script by Daniel Papp
% Wellcome Trust Centre for Neuroimaging (WTCN), London, UK.
% daniel.papp.13@ucl.ac.uk
% Adapted by Dr. Tobias Leutritz
% 
% Description:
% This script removes the coil sensitivty driven signal
% intensity modulation of MPM images. This is done as to correct for the
% different modulation of scans in case of inter-scan motion. 
% Also the case of only one acquistion per session is handled accordingly.
% 
% Outputs:
% corrected echoes, called sMT_manualnorm_echo-1, and the appropriate
% three "sensitivty maps", called MT_32ch_over_BC, or appropriate PD/T1
% 
% Reference: 
% D. Papp et al.: "Correction of Inter-Scan Motion Artifacts in
% Quantitative R1 Mapping by Accounting for Receive Coil Sensitivity
% Effects", MRM 2015 DOI 10.1002/mrm.26058

fprintf(1,'\n---------------- RF SENSITIVITY CORRECTION ----------------\n');

%==========================================================================
% Define processing parameters, defaults, input files...
%==========================================================================
rfsens_params = get_rfsens_params(jobsubj);

% for convenience:
calcpath = rfsens_params.calcpath;
supplpath = rfsens_params.supplpath;
smooth_kernel = rfsens_params.smooth_kernel;
json = rfsens_params.json;

% Input sensitivity maps: a pair of head coil/body coil sensitivity maps,
% either acquired once for the whole protocol, or once per contrast (MT,
% PD, T1):
MT_sensmaps = rfsens_params.input.MT.sensimages.fnames;
PD_sensmaps = rfsens_params.input.PD.sensimages.fnames;
T1_sensmaps = rfsens_params.input.T1.sensimages.fnames;

% Input MTw, PDw, T1w multiecho images:
MT_structurals = rfsens_params.input.MT.structurals.fnames;
PD_structurals = rfsens_params.input.PD.structurals.fnames;
T1_structurals = rfsens_params.input.T1.structurals.fnames;

%==========================================================================
% Coregistering the images: each sensmap onto the corresponding structural
%==========================================================================
MT_coregmaps = coreg_sens_to_struct_images(MT_structurals(1,:), MT_sensmaps, calcpath, 'MT');
PD_coregmaps = coreg_sens_to_struct_images(PD_structurals(1,:), PD_sensmaps, calcpath, 'PD');
T1_coregmaps = coreg_sens_to_struct_images(T1_structurals(1,:), T1_sensmaps, calcpath, 'T1');

%==========================================================================
% Smoothing the coregistered images
%==========================================================================
MT_smoothedmaps = smooth_sens_images(MT_coregmaps, smooth_kernel);
PD_smoothedmaps = smooth_sens_images(PD_coregmaps, smooth_kernel);
T1_smoothedmaps = smooth_sens_images(T1_coregmaps, smooth_kernel);


%==========================================================================
% Calculate quantitative RF sensitivity maps (HC/BC division)
%==========================================================================
MT_qsensmap = spm_imcalc(MT_smoothedmaps, fullfile(calcpath, 'sensMap_HC_over_BC_division_MT.nii'), 'i1./i2');MT_qsensmap = MT_qsensmap.fname;
PD_qsensmap = spm_imcalc(PD_smoothedmaps, fullfile(calcpath, 'sensMap_HC_over_BC_division_PD.nii'), 'i1./i2');PD_qsensmap = PD_qsensmap.fname;
T1_qsensmap = spm_imcalc(T1_smoothedmaps, fullfile(calcpath, 'sensMap_HC_over_BC_division_T1.nii'), 'i1./i2');T1_qsensmap = T1_qsensmap.fname;
% set metadata MT
input_files = MT_sensmaps;
Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
Output_hdr.history.output.imtype = 'Quantitative RF sensitivity map (HC/BC) for MTw images';
set_metadata(MT_qsensmap,Output_hdr,json);
% set metadata PD
input_files = PD_sensmaps;
Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
Output_hdr.history.output.imtype = 'Quantitative RF sensitivity map (HC/BC) for PDw images';
set_metadata(PD_qsensmap,Output_hdr,json);
% set metadata T1
input_files = T1_sensmaps;
Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
Output_hdr.history.output.imtype = 'Quantitative RF sensitivity map (HC/BC) for T1w images';
set_metadata(T1_qsensmap,Output_hdr,json);


%==========================================================================
% Normalise all input multi-echo images using the p.u. sensitivity maps
%==========================================================================
nMT = size(MT_structurals,1);
nPD = size(PD_structurals,1);
nT1 = size(T1_structurals,1);
MT_corrected_structurals = cell(nMT,1);
PD_corrected_structurals = cell(nPD,1);
T1_corrected_structurals = cell(nT1,1);

for i=1:nMT
    MT_corrected_structurals{i} = fullfile(calcpath, spm_file(spm_file(MT_structurals(i,:),'filename'),'suffix','_RFSC'));
    spm_imcalc({MT_structurals(i,:), MT_qsensmap}, MT_corrected_structurals{i}, 'i1./i2');
    %set metadata
    input_files = char(MT_structurals(i,:), MT_qsensmap);
    Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
    Output_hdr.history.output.imtype = 'RF sensitivity corrected MT-weighted echo';
    Output_hdr.history.output.units = 'a.u.';
    Output_hdr.acqpar = struct('RepetitionTime',get_metadata_val(MT_structurals(i,:),'RepetitionTime'), ...
        'EchoTime',get_metadata_val(MT_structurals(i,:),'EchoTime'), ...
        'FlipAngle',get_metadata_val(MT_structurals(i,:),'FlipAngle'));
    set_metadata(MT_corrected_structurals{i},Output_hdr,json);
end

for i=1:nPD
    PD_corrected_structurals{i} = fullfile(calcpath, spm_file(spm_file(PD_structurals(i,:),'filename'),'suffix','_RFSC'));
    spm_imcalc({PD_structurals(i,:), PD_qsensmap}, PD_corrected_structurals{i}, 'i1./i2');
    %set metadata
    input_files = char(PD_structurals(i,:), PD_qsensmap);
    Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
    Output_hdr.history.output.imtype = 'RF sensitivity corrected PD-weighted echo';
    Output_hdr.history.output.units = 'a.u.';
    Output_hdr.acqpar = struct('RepetitionTime',get_metadata_val(PD_structurals(i,:),'RepetitionTime'), ...
        'EchoTime',get_metadata_val(PD_structurals(i,:),'EchoTime'), ...
        'FlipAngle',get_metadata_val(PD_structurals(i,:),'FlipAngle'));
    set_metadata(PD_corrected_structurals{i},Output_hdr,json);
end

for i=1:nT1
    T1_corrected_structurals{i} = fullfile(calcpath, spm_file(spm_file(T1_structurals(i,:),'filename'),'suffix','_RFSC'));
    spm_imcalc({T1_structurals(i,:), T1_qsensmap}, T1_corrected_structurals{i}, 'i1./i2');
    %set metadata
    input_files = char(T1_structurals(i,:), T1_qsensmap);
    Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
    Output_hdr.history.output.imtype = 'RF sensitivity corrected T1-weighted echo';
    Output_hdr.history.output.units = 'a.u.';
    Output_hdr.acqpar = struct('RepetitionTime',get_metadata_val(T1_structurals(i,:),'RepetitionTime'), ...
        'EchoTime',get_metadata_val(T1_structurals(i,:),'EchoTime'), ...
        'FlipAngle',get_metadata_val(T1_structurals(i,:),'FlipAngle'));
    set_metadata(T1_corrected_structurals{i},Output_hdr,json);
end

%==========================================================================
% Finalise output
%==========================================================================
% copy quantitative RF sensitivity maps in Supplementary folder
copyfile(MT_qsensmap,fullfile(supplpath,spm_file(MT_qsensmap,'filename')));
try copyfile([spm_str_manip(MT_qsensmap,'r') '.json'],fullfile(supplpath,[spm_file(MT_qsensmap,'basename'), '.json'])); end %#ok<*TRYNC>
copyfile(PD_qsensmap,fullfile(supplpath,spm_file(PD_qsensmap,'filename')));
try copyfile([spm_str_manip(PD_qsensmap,'r') '.json'],fullfile(supplpath,[spm_file(PD_qsensmap,'basename'), '.json'])); end %#ok<*TRYNC>
copyfile(T1_qsensmap,fullfile(supplpath,spm_file(T1_qsensmap,'filename')));
try copyfile([spm_str_manip(T1_qsensmap,'r') '.json'],fullfile(supplpath,[spm_file(T1_qsensmap,'basename'), '.json'])); end %#ok<*TRYNC>

% save RF sensitivity processing parameters
spm_jsonwrite(fullfile(supplpath,'MPM_map_creation_rfsens_params.json'),rfsens_params,struct('indent','\t'));

% substitute the corrected maps to the output structure
jobsubj.raw_mpm.MT = char(MT_corrected_structurals);
jobsubj.raw_mpm.PD = char(PD_corrected_structurals);
jobsubj.raw_mpm.T1 = char(T1_corrected_structurals);

end


%% =======================================================================%
% Sort out all parameters required for the RFsens calculation.
%=========================================================================%
function rfsens_params = get_rfsens_params(jobsubj)

rfsens_params.json = hmri_get_defaults('json');
rfsens_params.calcpath = jobsubj.path.rfsenspath;
rfsens_params.respath = jobsubj.path.respath;
rfsens_params.supplpath = jobsubj.path.supplpath;
rfsens_params.smooth_kernel = hmri_get_defaults('RFsens.smooth_kernel');

% Input sensitivity maps: a pair of head coil/body coil sensitivity maps,
% either acquired once for the whole protocol, or once per contrast (MT,
% PD, T1). We first copy them to the calcpath directory to avoid working on
% raw data: 
if isfield(jobsubj.sensitivity,'RF_once')
    for i=1:length(jobsubj.sensitivity.RF_once)
        tmprawfnam = spm_file(jobsubj.sensitivity.RF_once{i},'number','');
        tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename')); %#ok<AGROW>
        copyfile(tmprawfnam, tmpfnam{i});
        try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
    end
    rfsens_params.input.MT.sensimages.fnames = char(tmpfnam);
    rfsens_params.input.PD.sensimages.fnames = char(tmpfnam);
    rfsens_params.input.T1.sensimages.fnames = char(tmpfnam);
    rfsens_params.senstype = 'RF_once: single sensitivity data set acquired for the whole MPM protocol';
elseif isfield(jobsubj.sensitivity,'RF_MPM')
    for i=1:length(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_MT)
        tmprawfnam = spm_file(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_MT{i},'number','');
        tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename')); %#ok<AGROW>
        copyfile(tmprawfnam, tmpfnam{i});
        try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
    end
    rfsens_params.input.MT.sensimages.fnames = char(tmpfnam);
    for i=1:length(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_PD)
        tmprawfnam = spm_file(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_PD{i},'number','');
        tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename'));
        copyfile(tmprawfnam, tmpfnam{i});
        try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
    end
    rfsens_params.input.PD.sensimages.fnames = char(tmpfnam);
    for i=1:length(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_T1)
        tmprawfnam = spm_file(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_T1{i},'number','');
        tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename')); 
        copyfile(tmprawfnam, tmpfnam{i});
        try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
    end
    rfsens_params.input.T1.sensimages.fnames = char(tmpfnam);  
    rfsens_params.senstype = 'RF_MPM: one sensitivity data set acquired per MPM contrast';
else
    error('RF sensitivity correction: no RF sensitivity data provided.');
end

% Input MTw, PDw, T1w multiecho images:
rfsens_params.input.MT.structurals.fnames = spm_file(char(jobsubj.raw_mpm.MT),'number','');
rfsens_params.input.PD.structurals.fnames = spm_file(char(jobsubj.raw_mpm.PD),'number','');
rfsens_params.input.T1.structurals.fnames = spm_file(char(jobsubj.raw_mpm.T1),'number','');

end


%% =======================================================================%
% To arrange the metadata structure for RFsens calculation output (can be
% adapted for each output separately, here are the most common defaults). 
%=========================================================================%
function metastruc = init_rfsens_output_metadata(input_files, rfsens_params)

proc.descrip = 'RF sensitivity correction';
proc.version = hmri_get_version;
proc.params = rfsens_params;
output.imtype = 'sensitivity map';
output.units = 'p.u.';
metastruc = init_output_metadata_structure(input_files, proc, output);

end

%% =======================================================================%
% Smoothing images 
%=========================================================================%
function smoosensfnam = smooth_sens_images(inputfnam, smookernel)
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = inputfnam;
matlabbatch{1}.spm.spatial.smooth.fwhm = smookernel*[1 1 1];
matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('smooth%d_', smookernel);
output_list = spm_jobman('run',matlabbatch);
smoosensfnam = output_list{1}.files;
clear matlabbatch
end

%% =======================================================================%
% Coregister sensitivity images onto structural images. Coregistration is
% done for HeadCoil and BodyCoil images separately. A tag corresponding to
% the contrast (MT/PD/T1) is added to avoid overwriting data when the same
% sensitivity map is used for each contrast.
%=========================================================================%
function coregsensfnam = coreg_sens_to_struct_images(strucfnam, sensfnam, calcpath, tag)
clear matlabbatch
for i = 1:size(sensfnam,1)
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {strucfnam(1,:)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sensfnam(i,:)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = sprintf('r%s_',tag);
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    output_list = spm_jobman('run',matlabbatch);
    coregsensfnam{i} = fullfile(calcpath, spm_file(output_list{1}.rfiles{1},'filename')); %#ok<AGROW>
end
coregsensfnam = coregsensfnam'; % must be a nx1 cellstr
clear matlabbatch
end