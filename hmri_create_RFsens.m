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

% Proceeds for each contrast
for ccon = 1:rfsens_params.ncon
    
    % Input sensitivity maps: a pair of head coil/body coil sensitivity maps,
    % either acquired once for the whole protocol, or once per contrast (MT,
    % PD, T1):
    sensmaps = rfsens_params.input(ccon).sensfnam;
    % Input MTw, PDw, T1w multiecho images:
    structurals = rfsens_params.input(ccon).fnam;
    
    %==========================================================================
    % Coregistering the images: each sensmap onto the corresponding structural
    %==========================================================================
    coregmaps = coreg_sens_to_struct_images(structurals(1,:), sensmaps, calcpath, rfsens_params.input(ccon).tag);
    
    %==========================================================================
    % Smoothing the coregistered images
    %==========================================================================
    smoothedmaps = smooth_sens_images(coregmaps, smooth_kernel);
    
    %==========================================================================
    % Calculate quantitative RF sensitivity maps (HC/BC division)
    %==========================================================================
    qsensmap = spm_imcalc(smoothedmaps, fullfile(calcpath, sprintf('sensMap_HC_over_BC_division_%s.nii',rfsens_params.input(ccon).tag)), 'i1./i2');
    qsensmap = qsensmap.fname;
    % set metadata
    input_files = sensmaps;
    Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
    Output_hdr.history.output.imtype = sprintf('Quantitative RF sensitivity map (HC/BC) for %sw images',rfsens_params.input(ccon).tag);
    set_metadata(qsensmap{con},Output_hdr,json);
    
    %==========================================================================
    % Normalise all input multi-echo images using the p.u. sensitivity maps
    %==========================================================================
    nSTRUCT = size(structurals,1);
    corrected_structurals = cell(nSTRUCT,1);
    
    for i=1:nSTRUCT
        corrected_structurals{i} = fullfile(calcpath, spm_file(spm_file(structurals(i,:),'filename'),'suffix','_RFSC'));
        spm_imcalc({structurals(i,:), qsensmap}, corrected_structurals{i}, 'i1./i2');
        % set metadata
        input_files = char(structurals(i,:), qsensmap);
        Output_hdr = init_rfsens_output_metadata(input_files, rfsens_params);
        Output_hdr.history.output.imtype = sprintf('RF sensitivity corrected %s-weighted echo',rfsens_params.input(ccon).tag);
        Output_hdr.history.output.units = 'a.u.';
        Output_hdr.acqpar = struct('RepetitionTime',get_metadata_val(structurals(i,:),'RepetitionTime'), ...
            'EchoTime',get_metadata_val(structurals(i,:),'EchoTime'), ...
            'FlipAngle',get_metadata_val(structurals(i,:),'FlipAngle'));
        set_metadata(corrected_structurals{i},Output_hdr,json);
    end
    
    %==========================================================================
    % Finalise output
    %==========================================================================
    % copy quantitative RF sensitivity maps in Supplementary folder
    copyfile(qsensmap,fullfile(supplpath,spm_file(qsensmap,'filename')));
    try copyfile([spm_str_manip(qsensmap,'r') '.json'],fullfile(supplpath,[spm_file(qsensmap,'basename'), '.json'])); end %#ok<*TRYNC>
    % substitute the corrected maps to the output structure
    jobsubj.raw_mpm.(rfsens_params.input(ccon).tag) = char(corrected_structurals);
    
end

% save RF sensitivity processing parameters
spm_jsonwrite(fullfile(supplpath,'MPM_map_creation_rfsens_params.json'),rfsens_params,struct('indent','\t'));

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

% Input structurals: determine which contrasts are available
ccon = 0;
% 1) try MTw contrast:
tmpfnam   = spm_file(char(jobsubj.raw_mpm.MT),'number',''); % P_mtw
if isempty(tmpfnam)
    rfsens_params.MTidx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    rfsens_params.input(ccon).fnam = tmpfnam;
    rfsens_params.input(ccon).tag = 'MT';
    rfsens_params.MTidx = ccon;
end
% 2) try PDw contrast:
tmpfnam   = spm_file(char(jobsubj.raw_mpm.PD),'number',''); % P_pdw
if isempty(tmpfnam)
    rfsens_params.PDidx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    rfsens_params.input(ccon).fnam = tmpfnam;
    rfsens_params.input(ccon).tag = 'PD';
    rfsens_params.PDidx = ccon;
end
% 3) try T1w contrast:
tmpfnam   = spm_file(char(jobsubj.raw_mpm.T1),'number',''); % P_t1w
if isempty(tmpfnam)
    rfsens_params.T1idx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    rfsens_params.input(ccon).fnam = tmpfnam;
    rfsens_params.input(ccon).tag = 'T1';
    rfsens_params.T1idx = ccon; % zero index means no contrast available
end
rfsens_params.ncon = ccon; % number of contrasts available
clear tmpfnam;

% Input sensitivity maps: a pair of head coil/body coil sensitivity maps,
% either acquired once for the whole protocol, or once per contrast (MT,
% PD, T1). We first copy them to the calcpath directory to avoid working on
% raw data:
if isfield(jobsubj.sensitivity,'RF_once')
    for i=1:length(jobsubj.sensitivity.RF_once)
        tmprawfnam = spm_file(jobsubj.sensitivity.RF_once{i},'number','');
        tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename'));
        copyfile(tmprawfnam, tmpfnam{i});
        try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
    end
    for ccon = 1:rfsens_params.ncon
        rfsens_params.input(ccon).sensfnam = char(tmpfnam);
    end
    rfsens_params.senstype = 'RF_once: single set of RF sensitivity maps acquired for all contrasts';
elseif isfield(jobsubj.sensitivity,'RF_per_contrast')
    if (rfsens_params.MTidx)
        for i=1:length(jobsubj.sensitivity.RF_per_contrast.raw_sens_MT)
            tmprawfnam = spm_file(jobsubj.sensitivity.RF_per_contrast.raw_sens_MT{i},'number','');
            tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename'));
            copyfile(tmprawfnam, tmpfnam{i});
            try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
        end
        rfsens_params.input.MT.sensfnam = char(tmpfnam);
    end
    if (rfsens_params.PDidx)
        for i=1:length(jobsubj.sensitivity.RF_per_contrast.raw_sens_PD)
            tmprawfnam = spm_file(jobsubj.sensitivity.RF_per_contrast.raw_sens_PD{i},'number','');
            tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename'));
            copyfile(tmprawfnam, tmpfnam{i});
            try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
        end
        rfsens_params.input.PD.sensfnam = char(tmpfnam);
    end
    if (rfsens_params.T1idx)
        for i=1:length(jobsubj.sensitivity.RF_per_contrast.raw_sens_T1)
            tmprawfnam = spm_file(jobsubj.sensitivity.RF_per_contrast.raw_sens_T1{i},'number','');
            tmpfnam{i} = fullfile(rfsens_params.calcpath,spm_file(tmprawfnam,'filename'));
            copyfile(tmprawfnam, tmpfnam{i});
            try copyfile([spm_str_manip(tmprawfnam,'r') '.json'],[spm_str_manip(tmpfnam{i},'r') '.json']); end
        end
        rfsens_params.input.T1.sensfnam = char(tmpfnam);
    end
    rfsens_params.senstype = 'RF_per_contrast: one sensitivity data set acquired per contrast (i.e. T1/PD/MT-weighted images)';
else
    error('RF sensitivity correction: no RF sensitivity data provided.');
end

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