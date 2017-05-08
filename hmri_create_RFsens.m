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

% retrieve some defaults
outdir = jobsubj.path.rfsenspath;
smooth_kernel = hmri_get_defaults('RFsens.smooth_kernel');
json = hmri_get_defaults('json');

%% calculate the RF sensitivity maps for each modality
% first assign sensitivity maps, save the original file names for history, 
if isfield(jobsubj.sensitivity,'RF_once')
    MT_sensmaps0 = char(jobsubj.sensitivity.RF_once);
    PD_sensmaps0 = char(jobsubj.sensitivity.RF_once);
    T1_sensmaps0 = char(jobsubj.sensitivity.RF_once);    
elseif isfield(jobsubj.sensitivity,'RF_MPM')
    MT_sensmaps0 = char(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_MT);
    PD_sensmaps0 = char(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_PD);
    T1_sensmaps0 = char(jobsubj.sensitivity.RF_MPM.raw_sens3.raw_sens_T1);
end

% prepare filenames, i.e. remove ',1' from spm_select
MT_sensmaps01 = spm_file(MT_sensmaps0(1,:),'number','');
MT_sensmaps02 = spm_file(MT_sensmaps0(2,:),'number','');
MT_sensmaps0 = char(MT_sensmaps01,MT_sensmaps02);
PD_sensmaps01 = spm_file(PD_sensmaps0(1,:),'number','');
PD_sensmaps02 = spm_file(PD_sensmaps0(2,:),'number','');
PD_sensmaps0 = char(PD_sensmaps01,PD_sensmaps02);
T1_sensmaps01 = spm_file(T1_sensmaps0(1,:),'number','');
T1_sensmaps02 = spm_file(T1_sensmaps0(2,:),'number','');
T1_sensmaps0 = char(T1_sensmaps01,T1_sensmaps02);

% afterwards copy and rename in order to have a set of maps for each modality
% especially in the case of only one measure

[~,~,ext] = fileparts(MT_sensmaps0(1,:));
MT_sensmap1 = strcat(outdir,filesep,'MT_sens_head',ext);
copyfile(deblank(MT_sensmaps0(1,:)),MT_sensmap1);
[~,~,ext] = fileparts(MT_sensmaps0(2,:));
MT_sensmap2 = strcat(outdir,filesep,'MT_sens_body',ext);
copyfile(deblank(MT_sensmaps0(2,:)),MT_sensmap2);
MT_sensmaps = char(MT_sensmap1,MT_sensmap2);
% set and write metadata
Output_hdr = struct('history',struct('procstep',[],'input',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'B1- map calculation: copy raw data';
for ctr = 1:size(MT_sensmaps0,1)
    Output_hdr.history.input.filename = MT_sensmaps0(ctr,:);
    input_hdr = get_metadata(MT_sensmaps0(ctr,:));
    if ~isempty(input_hdr{1})
        Output_hdr.history.input.history = input_hdr{1}.history;
    else
        Output_hdr.history.input.history = 'No history available.';
    end
    set_metadata(MT_sensmaps(ctr,:),Output_hdr,json);
end

[~,~,ext] = fileparts(PD_sensmaps0(1,:));
PD_sensmap1 = strcat(outdir,filesep,'PD_sens_head',ext);
copyfile(deblank(PD_sensmaps0(1,:)),PD_sensmap1);
[~,~,ext] = fileparts(PD_sensmaps0(2,:));
PD_sensmap2 = strcat(outdir,filesep,'PD_sens_body',ext);
copyfile(deblank(PD_sensmaps0(2,:)),PD_sensmap2);
PD_sensmaps = char(PD_sensmap1,PD_sensmap2);
% set and write metadata
Output_hdr = struct('history',struct('procstep',[],'input',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'B1- map calculation: copy raw data';
for ctr = 1:size(PD_sensmaps0,1)
    Output_hdr.history.input.filename = PD_sensmaps0(ctr,:);
    input_hdr = get_metadata(PD_sensmaps0(ctr,:));
    if ~isempty(input_hdr{1})
        Output_hdr.history.input.history = input_hdr{1}.history;
    else
        Output_hdr.history.input.history = 'No history available.';
    end
    set_metadata(PD_sensmaps(ctr,:),Output_hdr,json);
end

[~,~,ext] = fileparts(T1_sensmaps0(1,:));
T1_sensmap1 = strcat(outdir,filesep,'T1_sens_head',ext);
copyfile(deblank(T1_sensmaps0(1,:)),T1_sensmap1);
[~,~,ext] = fileparts(T1_sensmaps0(2,:));
T1_sensmap2 = strcat(outdir,filesep,'T1_sens_body',ext);
copyfile(deblank(T1_sensmaps0(2,:)),T1_sensmap2);
T1_sensmaps = char(T1_sensmap1,T1_sensmap2);
% set and write metadata
Output_hdr = struct('history',struct('procstep',[],'input',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'B1- map calculation: copy raw data';
for ctr = 1:size(T1_sensmaps0,1)
    Output_hdr.history.input.filename = T1_sensmaps0(ctr,:);
    input_hdr = get_metadata(T1_sensmaps0(ctr,:));
    if ~isempty(input_hdr{1})
        Output_hdr.history.input.history = input_hdr{1}.history;
    else
        Output_hdr.history.input.history = 'No history available.';
    end
    set_metadata(T1_sensmaps(ctr,:),Output_hdr,json);
end
    
% assign appropriate structurals
MT_structurals = char(jobsubj.raw_mpm.MT);
PD_structurals = char(jobsubj.raw_mpm.PD);
T1_structurals = char(jobsubj.raw_mpm.T1);

% create cell for the coregistered images by using the old filenames 
% before coregistration, and add the r_ marker to the filename
MT_coregmaps = cell(size(MT_sensmaps,1),1);

% coregistering the images. Probably could be done more elegantly
clear matlabbatch
for i = 1:2
    spm_jobman('initcfg')
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {deblank(MT_structurals(1,:))};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {deblank(MT_sensmaps(i,:))};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    spm_jobman('run',matlabbatch)
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map calculation: coregister with structural';
    Output_hdr.history.procstep.procpar = matlabbatch{1};
    Output_hdr.history.input{1}.filename = MT_structurals(1,:);
    Output_hdr.history.input{2}.filename = MT_sensmaps(i,:);
    for ctr = 1:2
        input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'coregistered raw sensitivity map';
    Output_hdr.history.output.units = 'a.u.';
    namelimits = strfind(MT_sensmaps(i,:),filesep);
    MT_coregmaps{i,:} = strcat(MT_sensmaps(i,1:namelimits(end)),'r_',(MT_sensmaps(i,namelimits(end)+1:end)));
    set_metadata(MT_coregmaps{i,:},Output_hdr,json);
end
clear matlabbatch
clear namelimits;

PD_coregmaps = cell(size(PD_sensmaps,1),1);
for i = 1:2
    spm_jobman('initcfg')
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {deblank(PD_structurals(1,:))};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {deblank(PD_sensmaps(i,:))};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    spm_jobman('run',matlabbatch)
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map calculation: coregister with structural';
    Output_hdr.history.procstep.procpar = matlabbatch{1};
    Output_hdr.history.input{1}.filename = PD_structurals(1,:);
    Output_hdr.history.input{2}.filename = PD_sensmaps(i,:);
    for ctr = 1:2
        input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'coregistered raw sensitivity map';
    Output_hdr.history.output.units = 'a.u.';
    namelimits = strfind(PD_sensmaps(i,:),filesep);
    PD_coregmaps{i,:} = strcat(PD_sensmaps(i,1:namelimits(end)),'r_',(PD_sensmaps(i,namelimits(end)+1:end)));
    set_metadata(PD_coregmaps{i,:},Output_hdr,json);
end
clear matlabbatch
clear namelimits;

T1_coregmaps = cell(size(T1_sensmaps,1),1);
for i = 1:2
    spm_jobman('initcfg')
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {deblank(T1_structurals(1,:))};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {deblank(T1_sensmaps(i,:))};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    spm_jobman('run',matlabbatch)
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map calculation: coregister with structural';
    Output_hdr.history.procstep.procpar = matlabbatch{1};
    Output_hdr.history.input{1}.filename = T1_structurals(1,:);
    Output_hdr.history.input{2}.filename = T1_sensmaps(i,:);
    for ctr = 1:2
        input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'coregistered raw sensitivity map';
    Output_hdr.history.output.units = 'a.u.';
    namelimits = strfind(T1_sensmaps(i,:),filesep);
    T1_coregmaps{i,:} = strcat(T1_sensmaps(i,1:namelimits(end)),'r_',(T1_sensmaps(i,namelimits(end)+1:end)));
    set_metadata(T1_coregmaps{i,:},Output_hdr,json);
end
clear matlabbatch
clear namelimits;

% create cell for the smoothed and coregistered files by using the old 
% filenames before coregistration, and add the smooth marker to it
MT_smoothedmaps = cell(size(MT_coregmaps,1),1);

% smooth coregistered images
% Use the coregistered images, and call spm_smooth for them
clear matlabbatch
for i = 1:2
    matlabbatch{1}.spm.spatial.smooth.data = {deblank(MT_coregmaps{i,:})};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [smooth_kernel smooth_kernel smooth_kernel];
    matlabbatch{1}.spm.spatial.smooth.prefix = strcat('smooth',num2str(smooth_kernel),'_');
    spm_jobman('run',matlabbatch)
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map calculation: smooth coregistered raw image';
    Output_hdr.history.procstep.procpar = matlabbatch{1};
    Output_hdr.history.input{1}.filename = MT_coregmaps{i,:};
    input_hdr = get_metadata(Output_hdr.history.input{1}.filename);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{1}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{1}.history = 'No history available.';
    end
    Output_hdr.history.output.imtype = 'smoothed raw sensitivity map';
    Output_hdr.history.output.units = 'a.u.';
    namelimits = strfind(MT_coregmaps{i,:},filesep);
    MT_smoothedmaps{i,:} = strcat((MT_coregmaps{i}(1:namelimits(end))),...
        strcat('smooth',num2str(smooth_kernel),'_'),...
        (MT_coregmaps{i}(namelimits(end)+1:end)));
    set_metadata(MT_smoothedmaps{i,:},Output_hdr,json);
end

PD_smoothedmaps = cell(size(PD_coregmaps,1),1);
clear matlabbatch
for i = 1:2
    matlabbatch{1}.spm.spatial.smooth.data = {deblank(PD_coregmaps{i,:})};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [smooth_kernel smooth_kernel smooth_kernel];
    matlabbatch{1}.spm.spatial.smooth.prefix = strcat('smooth',num2str(smooth_kernel),'_');
    spm_jobman('run',matlabbatch)
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map calculation: smooth coregistered raw image';
    Output_hdr.history.procstep.procpar = matlabbatch{1};
    Output_hdr.history.input{1}.filename = PD_coregmaps{i,:};
    input_hdr = get_metadata(Output_hdr.history.input{1}.filename);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{1}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{1}.history = 'No history available.';
    end
    Output_hdr.history.output.imtype = 'smoothed raw sensitivity map';
    Output_hdr.history.output.units = 'a.u.';
    namelimits = strfind(PD_coregmaps{i,:},filesep);
    PD_smoothedmaps{i,:} = strcat((PD_coregmaps{i}(1:namelimits(end))),...
        strcat('smooth',num2str(smooth_kernel),'_'),...
        (PD_coregmaps{i}(namelimits(end)+1:end)));
    set_metadata(PD_smoothedmaps{i,:},Output_hdr,json);
end

T1_smoothedmaps = cell(size(T1_coregmaps,1),1);
clear matlabbatch
for i = 1:2
    matlabbatch{1}.spm.spatial.smooth.data = {deblank(T1_coregmaps{i,:})};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [smooth_kernel smooth_kernel smooth_kernel];
    matlabbatch{1}.spm.spatial.smooth.prefix = strcat('smooth',num2str(smooth_kernel),'_');
    spm_jobman('run',matlabbatch)
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map calculation: smooth coregistered raw image';
    Output_hdr.history.procstep.procpar = matlabbatch{1};
    Output_hdr.history.input{1}.filename = T1_coregmaps{i,:};
    input_hdr = get_metadata(Output_hdr.history.input{1}.filename);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{1}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{1}.history = 'No history available.';
    end
    Output_hdr.history.output.imtype = 'smoothed raw sensitivity map';
    Output_hdr.history.output.units = 'a.u.';
    namelimits = strfind(T1_coregmaps{i,:},filesep);
    T1_smoothedmaps{i,:} = strcat((T1_coregmaps{i}(1:namelimits(end))),...
        strcat('smooth',num2str(smooth_kernel),'_'),...
        (T1_coregmaps{i}(namelimits(end)+1:end)));
    set_metadata(T1_smoothedmaps{i,:},Output_hdr,json);
end

% Divide the 32ch image with the BC image
[Filepath,~,~] = fileparts(MT_smoothedmaps{1});
clear matlabbatch;
matlabbatch{1}.spm.util.imcalc.input = {deblank(MT_smoothedmaps{1}),deblank(MT_smoothedmaps{2})}';
matlabbatch{1}.spm.util.imcalc.output = 'MT_32ch_over_BC.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {Filepath};
matlabbatch{1}.spm.util.imcalc.expression = 'i1./i2';
spm_jobman('run',matlabbatch)
% set and write metadata
Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'B1- map calculation: MT';
Output_hdr.history.procstep.procpar = matlabbatch{1};
Output_hdr.history.input{1}.filename = matlabbatch{1}.spm.util.imcalc.input{1};
Output_hdr.history.input{2}.filename = matlabbatch{1}.spm.util.imcalc.input{2};
for ctr = 1:2
    input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{ctr}.history = 'No history available.';
    end
end
Output_hdr.history.output.imtype = 'sensitivity map for MT';
Output_hdr.history.output.units = 'p.u.';
MT_map = strcat(Filepath,filesep,matlabbatch{1}.spm.util.imcalc.output);
set_metadata(MT_map(1,:),Output_hdr,json);


% generate manually normalized images
for i = 1:size(MT_structurals,1)
    FileName = [strcat(Filepath,filesep,'sMT_manualnorm_echo-',num2str(i)),'.nii'];
    placeholder = spm_imcalc({deblank(MT_structurals(i,:)), deblank(MT_map)}, FileName, 'i1./i2');
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[],'acqpar',...
        struct('RepetitionTime',[],'EchoTime',[],'FlipAngle',[])));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map application to MT';
    Output_hdr.history.procstep.procpar = 'i1./i2';
    Output_hdr.history.input{1}.filename = MT_structurals(i,:);
    Output_hdr.history.input{2}.filename = MT_map;
    Output_hdr.history.acqpar.RepetitionTime = ...
        get_metadata_val(MT_structurals(i,:),'RepetitionTime');
    Output_hdr.history.acqpar.EchoTime = ...
        get_metadata_val(MT_structurals(i,:),'EchoTime');
    Output_hdr.history.acqpar.FlipAngle = ...
        get_metadata_val(MT_structurals(i,:),'FlipAngle');
    for ctr = 1:2
        input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'RF sensitivity corrected MT';
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(FileName,Output_hdr,json);
end    

% change nifti headers
originalMT = cell(size(MT_structurals,1),1);
for i = 1:size(originalMT,1)
    originalMT{i,1} = nifti(deblank(MT_structurals(i,:)));
end

correctedMT = cell(size(MT_structurals,1),1);
correctedMTs = cell(size(MT_structurals,1),1);
for i = 1:size(correctedMT,1)
    FileName = [strcat(Filepath,filesep,'sMT_manualnorm_echo-',num2str(i)),'.nii'];
    correctedMT{i,1} = nifti(FileName);
    correctedMT{i,1}.descrip = originalMT{i,1}.descrip;
    create(correctedMT{i,1});
    correctedMTs{i,1} = FileName;
end

%% calculating the same for PD
% Divide the 32ch image with the BC image
[Filepath,~,~] = fileparts(PD_smoothedmaps{1});
clear matlabbatch;
matlabbatch{1}.spm.util.imcalc.input = {deblank(PD_smoothedmaps{1}),deblank(PD_smoothedmaps{2})}';
matlabbatch{1}.spm.util.imcalc.output = 'PD_32ch_over_BC.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {Filepath};
matlabbatch{1}.spm.util.imcalc.expression = 'i1./i2';
spm_jobman('run',matlabbatch)
% set and write metadata
Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'B1- map calculation: PD';
Output_hdr.history.procstep.procpar = matlabbatch{1};
Output_hdr.history.input{1}.filename = matlabbatch{1}.spm.util.imcalc.input{1};
Output_hdr.history.input{2}.filename = matlabbatch{1}.spm.util.imcalc.input{2};
for ctr = 1:2
    input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{ctr}.history = 'No history available.';
    end
end
Output_hdr.history.output.imtype = 'sensitivity map for PD';
Output_hdr.history.output.units = 'p.u.';
PD_map = strcat(Filepath,filesep,matlabbatch{1}.spm.util.imcalc.output);
set_metadata(PD_map(1,:),Output_hdr,json);

% generate manually normalized images
for i = 1:size(PD_structurals,1)
    FileName = [strcat(Filepath,filesep,'sPD_manualnorm_echo-',num2str(i)),'.nii'];
    placeholder = spm_imcalc({deblank(PD_structurals(i,:)), deblank(PD_map)}, FileName, 'i1./i2');
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[],'acqpar',...
        struct('RepetitionTime',[],'EchoTime',[],'FlipAngle',[])));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map application to PD';
    Output_hdr.history.procstep.procpar = 'i1./i2';
    Output_hdr.history.input{1}.filename = PD_structurals(i,:);
    Output_hdr.history.input{2}.filename = PD_map;
    Output_hdr.history.acqpar.RepetitionTime = ...
        get_metadata_val(PD_structurals(i,:),'RepetitionTime');
    Output_hdr.history.acqpar.EchoTime = ...
        get_metadata_val(PD_structurals(i,:),'EchoTime');
    Output_hdr.history.acqpar.FlipAngle = ...
        get_metadata_val(PD_structurals(i,:),'FlipAngle');
    for ctr = 1:2
        input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'RF sensitivity corrected PD';
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(FileName,Output_hdr,json);
end    

% change nifti headers
originalPD = cell(size(PD_structurals,1),1);
for i = 1:size(originalPD,1)
    originalPD{i,1} = nifti(deblank(PD_structurals(i,:)));
end

correctedPD = cell(size(PD_structurals,1),1);
correctedPDs = cell(size(PD_structurals,1),1);
for i = 1:size(correctedPD,1)
    FileName = [strcat(Filepath,filesep,'sPD_manualnorm_echo-',num2str(i)),'.nii'];
    correctedPD{i,1} = nifti(FileName);
    correctedPD{i,1}.descrip = originalPD{i,1}.descrip;
    create(correctedPD{i,1});
    correctedPDs{i,1} = FileName;
end

%% same thing for T1
% Divide the 32ch image with the BC image
[Filepath,~,~] = fileparts(T1_smoothedmaps{1});
clear matlabbatch;
matlabbatch{1}.spm.util.imcalc.input = {deblank(T1_smoothedmaps{1}),deblank(T1_smoothedmaps{2})}';
matlabbatch{1}.spm.util.imcalc.output = 'T1_32ch_over_BC.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {Filepath};
matlabbatch{1}.spm.util.imcalc.expression = 'i1./i2';
spm_jobman('run',matlabbatch)
% set and write metadata
Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'B1- map calculation: T1';
Output_hdr.history.procstep.procpar = matlabbatch{1};
Output_hdr.history.input{1}.filename = matlabbatch{1}.spm.util.imcalc.input{1};
Output_hdr.history.input{2}.filename = matlabbatch{1}.spm.util.imcalc.input{2};
for ctr = 1:2
    input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{ctr}.history = 'No history available.';
    end
end
Output_hdr.history.output.imtype = 'sensitivity map for T1';
Output_hdr.history.output.units = 'p.u.';
T1_map = strcat(Filepath,filesep,matlabbatch{1}.spm.util.imcalc.output);
set_metadata(T1_map(1,:),Output_hdr,json);

% generate manually normalized images
for i = 1:size(T1_structurals,1)
    FileName = [strcat(Filepath,filesep,'sT1_manualnorm_echo-',num2str(i)),'.nii'];
    placeholder = spm_imcalc({deblank(T1_structurals(i,:)), deblank(T1_map)}, FileName, 'i1./i2');
    % set and write metadata
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[],'acqpar',...
        struct('RepetitionTime',[],'EchoTime',[],'FlipAngle',[])));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'B1- map application to T1';
    Output_hdr.history.procstep.procpar = 'i1./i2';
    Output_hdr.history.input{1}.filename = T1_structurals(i,:);
    Output_hdr.history.input{2}.filename = T1_map;
    Output_hdr.history.acqpar.RepetitionTime = ...
        get_metadata_val(T1_structurals(i,:),'RepetitionTime');
    Output_hdr.history.acqpar.EchoTime = ...
        get_metadata_val(T1_structurals(i,:),'EchoTime');
    Output_hdr.history.acqpar.FlipAngle = ...
        get_metadata_val(T1_structurals(i,:),'FlipAngle');
    for ctr = 1:2
        input_hdr = get_metadata(Output_hdr.history.input{ctr}.filename);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = 'RF sensitivity corrected T1';
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(FileName,Output_hdr,json);
end    

% change nifti headers
originalT1 = cell(size(T1_structurals,1),1);
for i = 1:size(originalT1,1)
    originalT1{i,1} = nifti(deblank(T1_structurals(i,:)));
end

correctedT1 = cell(size(T1_structurals,1),1);
correctedT1s = cell(size(T1_structurals,1),1);
for i = 1:size(correctedT1,1)
    FileName = [strcat(Filepath,filesep,'sT1_manualnorm_echo-',num2str(i)),'.nii'];
    correctedT1{i,1} = nifti(FileName);
    correctedT1{i,1}.descrip = originalT1{i,1}.descrip;
    create(correctedT1{i,1});
    correctedT1s{i,1} = FileName;
end

% assign the corrected maps to the output structure
jobsubj.raw_mpm.MT = correctedMTs;
jobsubj.raw_mpm.PD = correctedPDs;
jobsubj.raw_mpm.T1 = correctedT1s;

end


%% =======================================================================%
% Sort out all parameters required for the RFsens calculation.
%=========================================================================%
function rfsens_params = get_rfsens_params(jobsubj)

rfsens_params.json = hmri_get_defaults('json');
rfsens_params.calcpath = jobsubj.path.rfsenspath;
rfsens_params.respath = jobsubj.path.respath;
rfsens_params.supplpath = jobsubj.path.supplpath;
rfsens_params.proc = hmri_get_defaults('RFsens');

end

%% =======================================================================%
% To arrange the metadata structure for RFsens calculation output.
%=========================================================================%
function metastruc = init_rfsens_output_metadata(input_files, rfsens_params)

proc.descrip = 'RF sensitivity calculation';
proc.version = hmri_get_version;
proc.params = rfsens_params;

% must be defined on the spot, default values here
output.imtype = 'sensitivity map';
output.units = 'p.u.';

metastruc = init_output_metadata_structure(input_files, proc, output);

end
