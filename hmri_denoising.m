function varargout = hmri_denoising(jobsubj)

% retrieve effective acquisition & processing parameters
denoising_params = get_denoising_params(jobsubj);
protocolfield = fieldnames(jobsubj.denoisingtype);
denoising_protocol = protocolfield{1};

% execute the chosen denoising method and define output
switch denoising_protocol
    case 'lcpca_denoise'
        [output_mag, output_phase] = hmri_calc_lcpcadenoise(denoising_params);
        varargout{1} = output_mag;
        varargout{2} = output_phase;
end
end

%=========================================================================%
% Write denoising parameters to denoising_params.
%=========================================================================%
function denoising_params = get_denoising_params(jobsubj)

% get denoising method
dntypename = fieldnames(jobsubj.denoisingtype);
denoising_protocol = dntypename{1};

% init defaults filename and optional-defaults bool
deffnam = '';
custom_def = false;

% load customized defaults parameters from customized denoising defaults file if any
% (the customized defaults file must be run to overwrite the standard
% defaults parameters)
if isfield(jobsubj.denoisingtype.(denoising_protocol),'DNparameters')
    % first reinitialise processing parameters to standard defaults:
    hmri_denoising_defaults;
    deffnam = fullfile(fileparts(mfilename('fullpath')),'config','hmri_denoising_defaults.m');
    custom_def = false;

    % then, if customized defaults file available, run it to overwrite
    % standard defaults parameters.
    if isfield(jobsubj.denoisingtype.(denoising_protocol).DNparameters,'DNdefaults')
        deffnam = jobsubj.denoisingtype.(denoising_protocol).DNparameters.DNdefaults;
        spm('Run',deffnam);
        custom_def = true;
    end
end

% set all denoising defaults and set the defaults file to be true
denoising_params = hmri_get_defaults(['denoising.' denoising_protocol]);
denoising_params.defaults_file = deffnam;
denoising_params.custom_defaults = custom_def;

% flags for logging information and warnings
denoising_params.defflags = jobsubj.log.flags; % default flags
denoising_params.nopuflags = jobsubj.log.flags; % force no Pop-Up
denoising_params.nopuflags.PopUp = false;

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
denoising_params.SPMver = sprintf('%s (%s)', v, r);

% Load input data
denoising_params.mag_img = cellstr(char(spm_file(jobsubj.mag_img,'number','')));
denoising_params.phase_img = cellstr(char(spm_file(jobsubj.phase_img,'number','')));
denoising_params.phase_bool = any(~cellfun(@isempty, denoising_params.phase_img));
denoising_params.mag_bool = any(~cellfun(@isempty, denoising_params.mag_img));

% processing can continue if only magnitude images were entered but
% only warn that optional phase img are missing
if ~denoising_params.phase_bool || ~isfield(jobsubj,'phase_img')
    hmri_log('Warning: No (optional) phase images were entered, denoising will continue with only magnitude images', denoising_params.defflags);
end

% Denoising method-specific parameters
% Load all the batch entered and possibly user-modified parameters
switch denoising_protocol
    case 'lcpca_denoise'

        dnstruct = jobsubj.denoisingtype.lcpca_denoise;
        denoising_params.ngbsize = dnstruct.ngbsize;
        denoising_params.std = dnstruct.std;
        denoising_params.output_path = jobsubj.path.dnrespath;
        denoising_params.supp_path = jobsubj.path.supplpath;

        % Print lcpca denoising parameters ngbsize and std cut off
        print_lcpca_params.Neighborhood_Size = denoising_params.ngbsize;
        print_lcpca_params.Standard_Deviation_Cutoff = denoising_params.std;
        printdnstruct = printstruct(print_lcpca_params);
        hmri_log(sprintf('Lcpca Denoising Parameters:\n\n%s', ...
            printdnstruct),denoising_params.defflags);

end

end

%===============================================================================================%
% Calculate LCPCA-denoising: Java-Matlab interface for Lcpca-denoising
%=================================================================================================%
function [output_mag, output_phase] = hmri_calc_lcpcadenoise(lcpcadenoiseparams)

% Reference for Lcpca denoising (Java module originally written by Pilou Bazin):
%   Bazin, et al. (2019) "Denoising High-Field Multi-Dimensional MRI With Local
%   Complex PCA", Front. Neurosci. doi:10.3389/fnins.2019.01066

% get the path to the jar files
mfilepath = fileparts(mfilename("fullpath"));
jarcommons = fullfile(fullfile(mfilepath,'jar'),'commons-math3-3.5.jar');
jarmipav = fullfile(fullfile(mfilepath,'jar'),'Jama-mipav.jar');
jarlcpca = fullfile(fullfile(mfilepath,'jar'),'lcpca2.jar');

% add required .jar files to the path
javaaddpath(jarcommons)
javaaddpath(jarmipav)
javaaddpath(jarlcpca)

% Read processing parameters from input
image_list = cellstr(lcpcadenoiseparams.mag_img);
phase_list = cellstr(lcpcadenoiseparams.phase_img);
ngb_size = lcpcadenoiseparams.ngbsize;
stdev_cutoff = lcpcadenoiseparams.std;
min_dimension = lcpcadenoiseparams.min_dimension;
max_dimension = lcpcadenoiseparams.max_dimension;
unwrap = lcpcadenoiseparams.unwrap;
rescale_phs = lcpcadenoiseparams.rescale_phs;
process_2d = lcpcadenoiseparams.process_2d;
use_rmt = lcpcadenoiseparams.use_rmt;
output_path = cellstr(lcpcadenoiseparams.output_path);
supp_path = cellstr(lcpcadenoiseparams.supp_path);
lcpcaflags = lcpcadenoiseparams.defflags;

% Get the full input file list
phase_bool = ~any(cellfun(@isempty,phase_list));
if phase_bool
    fullim_list = [image_list; phase_list];
else
    fullim_list = image_list;
end

% set the metadata mod
json = hmri_get_defaults('json');

% Create denoising object
noiseObj= javaObject('nl.uva.imcn.algorithms.LocalComplexPCADenoising');

% Get the image_number and dimensions from the first image of magnitude images
imdatavol = spm_vol(image_list{1});
resolutions = sqrt(sum(imdatavol.mat(1:3,1:3).^2));
imdatamatrix = spm_read_vols(imdatavol);
dimensions = size(imdatamatrix);
image_number = length(image_list);

% Set the image_number, dimensions and resolution for processing
noiseObj.setImageNumber(image_number); % number of echoes

% Set dimensions based on whether it is 3D or 4D data
% Throw an error if neither
if length(dimensions) == 3
    noiseObj.setDimensions(dimensions(1),dimensions(2),dimensions(3));
elseif length(dimensions) == 4
    noiseObj.setDimensions(dimensions(1),dimensions(2),dimensions(3),dimensions(4));
else
    msg='The raw data does not have the correct dimensions (must be 3D or 4D data) for processing: please check your data!';
    hmri_log(msg, lcpcaflags);
    error(msg)
end
noiseObj.setResolutions(resolutions(1),resolutions(2),resolutions(3));

% Loop through all reshaped echos
for echo = 1:length(image_list)
    % read vol from filepath
    datavol = spm_vol(image_list{echo});
    datamatrix = spm_read_vols(datavol);
    % place vol for denoising
    noiseObj.setMagnitudeImageAt(echo-1, reshape(datamatrix, [1 numel(datamatrix)]));

    % Add phase images if they exist
    if phase_bool
        % read vol from filepath
        datavol = spm_vol(phase_list{echo});
        datamatrix = spm_read_vols(datavol);
        % place vol for denoising
        noiseObj.setPhaseImageAt(echo-1, reshape(datamatrix, [1 numel(datamatrix)]));
    end
end

% set all other params
noiseObj.setPatchSize(ngb_size)
noiseObj.setStdevCutoff(stdev_cutoff)
noiseObj.setMinimumDimension(min_dimension)
noiseObj.setMaximumDimension(max_dimension)
noiseObj.setUnwrapPhase(unwrap)
noiseObj.setRescalePhase(rescale_phs)
noiseObj.setProcessSlabIn2D(process_2d)
noiseObj.setRandomMatrixTheory(use_rmt)

% Execute with all the input parameters
% Catch errors thrown by the java script
lcpcaflags_nopopup = lcpcaflags;
lcpcaflags_nopopup.PopUp = false;
hmri_log(sprintf('Executing Lcpca-denoising (Java) module \n'), lcpcaflags_nopopup);
try
    noiseObj.execute;
catch ME
    disp(ME.message)
    if isa(ME,'matlab.exception.JavaException')
        hmri_log('There was an error in the Java code: could not execute Lcpca denoising!', lcpcaflags_nopopup);
        if strcmp(ME.ExceptionObject,'java.lang.OutOfMemoryError: Java heap space')
            msg = sprintf(['Error: Lcpca denoising did not run because the Java heap space is too small for the amount of data to be denoised.\n' ...
                           '       Please try increasing the java heap size in the Matlab preferences.\n']);
            hmri_log(msg, lcpcaflags_nopopup);
        end        
    end
    rethrow(ME)
end
hmri_log('Lcpca-denoising (Java) module executed successfully', lcpcaflags_nopopup);

% initialize the cells to be populated with fullpaths of output magnitude and phase images
out_mag = cell(1,length(image_list));
idx_mag = 1;
if phase_bool 
    out_phase = cell(1,length(image_list));
else
    out_phase = {};
end
idx_phase = 1;

% Get the results for all echos and reshape
for echo = 1:length(image_list)
    datamatrix = noiseObj.getDenoisedMagnitudeImageAt(echo-1);
    volumedata = reshape(datamatrix, dimensions);

    % Write the volume to .nii with an update to standard .nii header
    firstfile = image_list{echo};
    filehdr = spm_vol(image_list{echo});

    [path,filename,~] = fileparts(firstfile);
    [~,mainfilename,~] = fileparts(firstfile);
    filename = strcat('LcpcaDenoised_',filename,'.nii');

    outfname = fullfile(output_path{1}, filename);
    filehdr.fname = outfname;
    filehdr.descrip = strcat(filehdr.descrip, ' + lcpca denoised');
    spm_write_vol(filehdr, volumedata);

    % write metadata as extended header and sidecar json
    Output_hdr = init_dn_output_metadata(fullim_list, lcpcadenoiseparams);
    Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (LCPCA)'];
    Output_hdr.history.output.imtype = 'Denoised (LCPCA)';
    % add acquisition data if available (otherwise fields will be empty)
    jsonfilename = fullfile(path,strcat(mainfilename,'.json'));
    if exist(jsonfilename, 'file') ==2
        try
            jsondata = spm_jsonread(jsonfilename);
            data_RepetitionTime = get_metadata_val(jsondata,'RepetitionTime');
            data_EchoTime = get_metadata_val(jsondata,'EchoTime' );
            data_FlipAngle = get_metadata_val(jsondata, 'FlipAngle');
            Output_hdr.acqpar = struct('RepetitionTime',data_RepetitionTime, ...
                'EchoTime',data_EchoTime,'FlipAngle',data_FlipAngle);
        catch
            hmri_log('Although json sidecar file was found, the writing of acquisition metadata failed', lcpcaflags_nopopup);

        end
    else
        hmri_log('No json sidecar file was found, skipping the writing of acquisition metadata', lcpcaflags_nopopup);
    end

    % Set all the metadata
    set_metadata(outfname,Output_hdr,json);

    % add image to the output list
    out_mag{idx_mag} = outfname;
    idx_mag = idx_mag + 1;

    % Also write the phase images if applicable
    if phase_bool
        datamatrix = noiseObj.getDenoisedPhaseImageAt(echo-1);
        volumedata = reshape(datamatrix, dimensions);

        % Write the volume to .nii with an update to standard .nii header
        firstfile = phase_list{echo};
        filehdr = spm_vol(phase_list{echo});

        [path,filename,~] = fileparts(firstfile);
        [~,mainfilename,~] = fileparts(firstfile);
        filename = strcat('LcpcaDenoised_',filename,'.nii');

        outfname = fullfile(output_path{1}, filename);
        filehdr.fname = outfname;
        filehdr.descrip = strcat(filehdr.descrip, ' + lcpca denoised');
        spm_write_vol(filehdr, volumedata);

        % write metadata as extended header and sidecar json
        Output_hdr = init_dn_output_metadata(fullim_list, lcpcadenoiseparams);
        Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (LCPCA)'];
        Output_hdr.history.output.imtype = 'Denoised (LCPCA)';
        % add acquisition data if available (otherwise fields will be empty)
        jsonfilename = fullfile(path,strcat(mainfilename,'.json'));
        if exist(jsonfilename, 'file') ==2
            try
                jsondata = spm_jsonread(jsonfilename);
                data_RepetitionTime = get_metadata_val(jsondata,'RepetitionTime');
                data_EchoTime = get_metadata_val(jsondata,'EchoTime' );
                data_FlipAngle = get_metadata_val(jsondata, 'FlipAngle');
                Output_hdr.acqpar = struct('RepetitionTime',data_RepetitionTime, ...
                    'EchoTime',data_EchoTime,'FlipAngle',data_FlipAngle);
            catch
                hmri_log('Although json sidecar file was found, the writing of acquisition metadata failed', lcpcaflags_nopopup);

            end
        else
            hmri_log('No json sidecar file was found, skipping the writing of acquisition metadata', lcpcaflags_nopopup);
        end

        % Set all the metadata
        set_metadata(outfname,Output_hdr,json);

        % add image to the output list
        out_phase{idx_phase} = outfname;
        idx_phase = idx_phase + 1;

    end

end

% Take out computed outputs
output_mag = out_mag;
output_phase = out_phase;

% save estimated local dimensions and residuals (between input and denoised images)
dim_img = reshape(noiseObj.getLocalDimensionImage(), dimensions);
save(fullfile(supp_path{1}, 'dim_img.mat'), 'dim_img')
err_img = reshape(noiseObj.getNoiseFitImage(), dimensions);
save(fullfile(supp_path{1}, 'err_img.mat'), 'err_img')

% Clear object and remove .jar from path properly
clear("noiseObj")
javarmpath(jarcommons)
javarmpath(jarmipav)
javarmpath(jarlcpca)

end

%=========================================================================%
% To arrange the metadata structure for denoising output.
%=========================================================================%
function metastruc = init_dn_output_metadata(input_files, denoising_params)

proc.descrip = ['hMRI toolbox - ' mfilename '.m - Denoising'];
proc.version = hmri_get_version;
proc.params = denoising_params;

output.imtype = 'Denoised Image';
output.units = '';

metastruc = init_output_metadata_structure(input_files, proc, output);

end

%=========================================================================%
% To print a structure into text - assumes simple structure (no
% sub-structure in it at this point).
%=========================================================================%
function s = printstruct(struc)

s = '';
fntmp = fieldnames(struc);
for cf = 1:length(fntmp)
    s = sprintf('%s %16s: %s\n', s, fntmp{cf}, num2str(struc.(fntmp{cf})));
end
end
