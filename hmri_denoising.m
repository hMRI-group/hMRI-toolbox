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
    case 'mppca_denoise'
        [output_mag, output_phase] = hmri_calc_mppcadenoise(denoising_params);
        varargout{1} = output_mag;
        varargout{2} = output_phase;
    case 'paws_denoise'
        [output_mag, output_phase] = hmri_calc_adaptivedenoising(denoising_params);
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

% Load contrast tags
denoising_params.mag_img_cotags= string(jobsubj.mag_img_cotags);
denoising_params.phase_img_cotags= string(jobsubj.phase_img_cotags);

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
        case 'mppca_denoise'

        denoising_params.mask = jobsubj.denoisingtype.(denoising_protocol).mask;
        denoising_params.ngbsize = jobsubj.denoisingtype.(denoising_protocol).ngbsize;

        denoising_params.output_path = jobsubj.path.dnrespath;
        denoising_params.supp_path = jobsubj.path.supplpath;
        % Print mppca denoising parameters 
        print_mppca_params.ngbsize = denoising_params.ngbsize;
        printdnstruct = printstruct(print_mppca_params);
        hmri_log(sprintf('MP-PCA Denoising Parameters:\n\n%s', ...
            printdnstruct),denoising_params.defflags);

        case 'paws_denoise'
            
        denoising_params.mask = jobsubj.denoisingtype.(denoising_protocol).mask;
        denoising_params.kstar = jobsubj.denoisingtype.(denoising_protocol).kstar;
        denoising_params.patchsize = jobsubj.denoisingtype.(denoising_protocol).patchsize;
        denoising_params.ladjust = jobsubj.denoisingtype.(denoising_protocol).ladjust;
        denoising_params.output_path = jobsubj.path.dnrespath;
        denoising_params.supp_path = jobsubj.path.supplpath;

        % Get separate contrasts
        contrasts = {"pdw","t1w","mtw"};
        cotags= denoising_params.mag_img_cotags;
        full_mag_img=denoising_params.mag_img;
        for c = 1:length(contrasts)
            con = contrasts{c};
            con_filter = (cotags==con);
            denoising_params.(con) = full_mag_img(con_filter);
        end

        % Get TE for further processing
        denoising_params.t1w_te = get_te(char(denoising_params.t1w));
        denoising_params.pdw_te = get_te(char(denoising_params.pdw));
        denoising_params.mtw_te = get_te(char(denoising_params.mtw));

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
    dt = filehdr.dt;
    updated_descrip = strcat(filehdr.descrip, ' + lcpca denoised');
    denoised_nifti = hmri_create_nifti(outfname, filehdr, dt, updated_descrip);
    denoised_nifti.dat(:,:,:) = volumedata;

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
        dt = filehdr.dt;
        updated_descrip = strcat(filehdr.descrip, ' + lcpca denoised');
        denoised_nifti = hmri_create_nifti(outfname, filehdr, dt, updated_descrip);
        denoised_nifti.dat(:,:,:)=volumedata;

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

%%save estimated local dimensions and residuals (between input and denoised images)
%get metadata from first echo of magnitude images and write local dims as
%nifti
dim_img = reshape(noiseObj.getLocalDimensionImage(), dimensions);
filehdr = spm_vol(image_list{1});
filename = strcat('dim_img','.nii');
outfname = fullfile(supp_path{1}, filename);
filehdr.fname = outfname;
filehdr.descrip = 'local dimension map';
dt=[spm_type('float32'),spm_platform('bigend')];
Ni_dim = hmri_create_nifti(outfname, filehdr, dt, filehdr.descrip);
Ni_dim.dat(:,:,:)=dim_img;


%get metadata from first echo of magnitude images and write residual maps
%as nifti
err_img = reshape(noiseObj.getNoiseFitImage(), dimensions);
filehdr = spm_vol(image_list{1});
filename = strcat('err_img','.nii');
outfname = fullfile(supp_path{1}, filename);
filehdr.fname = outfname;
filehdr.descrip = 'residual maps';
dt=[spm_type('float32'),spm_platform('bigend')];
Ni_err = hmri_create_nifti(outfname, filehdr, dt, filehdr.descrip);
Ni_err.dat(:,:,:)=err_img;


% Clear object and remove .jar from path properly
clear("noiseObj")
javarmpath(jarcommons)
javarmpath(jarmipav)
javarmpath(jarlcpca)

end

%===============================================================================================%
% Calculate MP-PCA-denoising
%=================================================================================================%

% MP-PCA-Denoising paper references:
% Veraart, Jelle, et al. "Denoising of diffusion MRI using random matrix theory." Neuroimage 142 (2016): 394-406.
% Does, Mark D., et al. "Evaluation of principal component analysis image denoising on multi‐exponential MRI relaxometry." Magnetic resonance in medicine 81.6 (2019): 3503-3514.


function [output_mag, output_phase] = hmri_calc_mppcadenoise(mppcadenoiseparams)

% define the flag for hmri log
mppcaflags = mppcadenoiseparams.defflags;
mppcaflags_nopopup = mppcaflags;
mppcaflags_nopopup.PopUp = false;

% set the metadata mod
json = hmri_get_defaults('json');

% read from the input the processing parameters
image_list = cellstr(mppcadenoiseparams.mag_img);
firstIm= image_list{1};

% get params from 1st image and init variables
image = spm_vol(firstIm);
imagevol = spm_read_vols(image);
imsize = size(imagevol);
imlen = length(image_list);
phase_list = cellstr(mppcadenoiseparams.phase_img);
phscale=1;

% if phase images are entered prepare images for further processing
if ~isempty(phase_list{1})
% takes magnitude and phase images in
% bundles them together realPart+complexPart into array
% feeds them to mppca-denoising to get get output mag, phase images

% init image list for further processing
imglist = zeros(imsize(1), imsize(2), imsize(3), imlen);

    for i=1:length(image_list)
        mag_imgstr = spm_vol(image_list{i});
        mag_imgvol = spm_read_vols(mag_imgstr);

        phase_imgstr = spm_vol(phase_list{i});
        phase_imgvol = spm_read_vols(phase_imgstr);

        complex_vol = mag_imgvol.*(exp((1i).*(phscale*phase_imgvol)));
        imglist(:,:,:,i)=complex_vol;
    end
    fulldatamat = imglist;
    image_list=cat(1,image_list,phase_list);
else
    % process and reformat images for MPPCA
    fulldatamat = zeros(imsize(1), imsize(2), imsize(3), imlen);
    for ii = 1:imlen
    currentIm = image_list{ii};    
    image = spm_vol(currentIm);
    imagevol = spm_read_vols(image);
    fulldatamat(:,:,:,ii) = imagevol;
    end
end

% get and set mppca denoising params
ngb_size = mppcadenoiseparams.ngbsize;
mask = mppcadenoiseparams.mask;
if (iscell(mask) && isempty(mask{1})) || isempty(mask)
    mask = [];
else
    mask_str = spm_vol(mask);
    mask_vols = spm_read_vols(mask_str{1});
    dims= size(mask_vols);
    premask = true([dims 1]);
    premask(:,:,:,1)= mask_vols;
    premask = logical(premask);
    mask = premask;
end
window = [ngb_size ngb_size ngb_size];
output_path = cellstr(mppcadenoiseparams.output_path);
supp_path = cellstr(mppcadenoiseparams.supp_path);

% Execute MPPCA 
hmri_log(sprintf('Executing MPPCA denoising \n'), mppcaflags_nopopup);
try
    % apply mppca denoising take out and set variables
[dn_image, map_noise_var, num_pc] = mppca_denoise(fulldatamat, window, mask);
catch ME
    % write the error description to hmrilog and rethrow to matlab
    err_str=['There was an internal error in the MPPCA code: ' ME.message];
    hmri_log(err_str, mppcaflags_nopopup);
    rethrow(ME)
end
hmri_log('MPPCA denoising executed successfully', mppcaflags_nopopup);

% recalculate magnitude and phase images in case of phase input
if ~isempty(phase_list{1})
    magimg = abs(dn_image);
    phimg = angle(dn_image);
else
    magimg= dn_image;
end


% init index
idx_mag=1;

% get the results for all echos and reshape
for echo = 1:imlen
    % get denoised volume
    volumedata  = magimg(:,:,:,echo);
    % write the volume to .nii with an update to standard .nii header
    firstfile = image_list{echo};
    filehdr = spm_vol(image_list{echo});

    [path,filename,ext] = fileparts(firstfile);
    [~,mainfilename,~] = fileparts(firstfile);
    filename = strcat('MppcaDenoised_',filename,'.nii');

    outfname = fullfile(output_path{1}, filename);
    dt = filehdr.dt;
    updated_descrip = strcat(filehdr.descrip, ' + mppca denoised');
    denoised_nifti = hmri_create_nifti(outfname, filehdr, dt, updated_descrip);
    denoised_nifti.dat(:,:,:) = volumedata;

    % write metadata as extended header and sidecar json
    Output_hdr = init_dn_output_metadata(image_list, mppcadenoiseparams);
    Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (MP-PCA)'];
    Output_hdr.history.output.imtype = 'Denoising (MP-PCA)';
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
            hmri_log('Although json sidecar file were found, the writing of acquisition metadata failed', mppcaflags_nopopup);  
      
       end
    else
      hmri_log('No json sidecar file were found, skipping the writing of acquisition metadata', mppcaflags_nopopup);  
    end

    % set all the metadata
    set_metadata(outfname,Output_hdr,json);

    % add image to the output list
    out_mag{idx_mag} = outfname;
    idx_mag = idx_mag +1;


end
% set output
output_mag = out_mag;

output_phase= [];
% process further if phase images are entered
if ~isempty(phase_list{1})
idx_phase=1;
% get the results for all echos and reshape
for echo = 1:imlen
    % get denoised volume
    volumedata  = phimg(:,:,:,echo);
    % write the volume to .nii with an update to standard .nii header
    firstfile = phase_list{echo};
    filehdr = spm_vol(phase_list{echo});

    [path,filename,ext] = fileparts(firstfile);
    [~,mainfilename,~] = fileparts(firstfile);
    filename = strcat('MppcaDenoised_',filename,'.nii');

    outfname = fullfile(output_path{1}, filename);
    dt = filehdr.dt;
    updated_descrip = strcat(filehdr.descrip, ' + mppca denoised');
    denoised_nifti = hmri_create_nifti(outfname, filehdr, dt, updated_descrip);
    denoised_nifti.dat(:,:,:)=volumedata;

    % write metadata as extended header and sidecar json
    Output_hdr = init_dn_output_metadata(phase_list, mppcadenoiseparams);
    Output_hdr.history.procstep.descrip = [Output_hdr.history.procstep.descrip ' (MP-PCA)'];
    Output_hdr.history.output.imtype = 'Denoising (MP-PCA)';
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
            hmri_log('Although json sidecar file were found, the writing of acquisition metadata failed', mppcaflags_nopopup);  
      
       end
    else
      hmri_log('No json sidecar file were found, skipping the writing of acquisition metadata', mppcaflags_nopopup);  
    end

    % set all the metadata
    set_metadata(outfname,Output_hdr,json);

    % add image to the output list
    out_phase{idx_phase} = outfname;
    idx_phase = idx_phase +1;


end

% set output
output_phase =out_phase;
end

% save map_noise_var (map of estimated noise variance) as a nifti image to Results/Supplementary
filehdr = spm_vol(firstIm);
filename = strcat('map_noise_var','.nii');
outfname = fullfile(supp_path{1}, filename);
filehdr.fname = outfname;
filehdr.descrip = 'map of estimated noise variance';
dt=[spm_type('float32'),spm_platform('bigend')];
Ni_dim = hmri_create_nifti(outfname, filehdr, dt, filehdr.descrip);
Ni_dim.dat(:,:,:)=map_noise_var;

% save num_pc (number of detected signal principal components) as a nifti
% image to Results/Supplementary
filehdr = spm_vol(firstIm);
filename = strcat('num_pc','.nii');
outfname = fullfile(supp_path{1}, filename);
filehdr.fname = outfname;
filehdr.descrip = 'number of detected signal principal components';
dt=[spm_type('float32'),spm_platform('bigend')];
Ni_dim = hmri_create_nifti(outfname, filehdr, dt, filehdr.descrip);
Ni_dim.dat(:,:,:)=num_pc;

end


%===============================================================================================%
% Calculate adaptive-denoising
%=================================================================================================%

function [output_mag, output_phase] = hmri_calc_adaptivedenoising(denoising_params)

%define the flag for log
pawsflags = denoising_params.defflags;
pawsflags_nopopup = pawsflags;
pawsflags_nopopup.PopUp = false;

%%Read from the input the processing parameters
image_list = cellstr(denoising_params.mag_img);
firstIm= image_list{1};
%Get params from 1st image and init variables
image = spm_vol(firstIm);
imagevol = spm_read_vols(image);
imsize = size(imagevol);
imlen = length(image_list);
wghts = sqrt(sum(image.mat(1:3,1:3).^2));

%%Read pawsdenoise params
mask = denoising_params.mask;
kstar = denoising_params.kstar;
patchsize = denoising_params.patchsize;
ladjust = denoising_params.ladjust;

alpha= denoising_params.alpha;
fit_method = denoising_params.fit_method;
mscbw = denoising_params.mscbw;
spmin = denoising_params.spmin;
lambda0 = denoising_params.lambda0;
corr_fac_patchsize = denoising_params.corr_fac_patchsize;

% prepare weighted_data for further processing

weighted_data(1).data = zeros(imsize(1),imsize(2),imsize(3),length(denoising_params.t1w));
TE_array = [denoising_params.t1w_te.te];
weighted_data(1).TE = TE_array;


for img=1:length(denoising_params.t1w)
    im_str = spm_vol(denoising_params.t1w(img));
    im_vol = spm_read_vols(im_str{1}); 
    weighted_data(1).data(:,:,:,img) = im_vol;

end

weighted_data(2).data = zeros(imsize(1),imsize(2),imsize(3),length(denoising_params.mtw));
TE_array = [denoising_params.mtw_te.te];
weighted_data(2).TE = TE_array;

for img=1:length(denoising_params.mtw)
    im_str = spm_vol(denoising_params.mtw(img));
    im_vol = spm_read_vols(im_str{1}); 
    weighted_data(2).data(:,:,:,img) = im_vol;

end

weighted_data(3).data = zeros(imsize(1),imsize(2),imsize(3),length(denoising_params.pdw));
TE_array = [denoising_params.pdw_te.te];
weighted_data(3).TE = TE_array;

for img=1:length(denoising_params.pdw)
    im_str = spm_vol(denoising_params.pdw(img));
    im_vol = spm_read_vols(im_str{1}); 
    weighted_data(3).data(:,:,:,img) = im_vol;
end

end

%=========================================================================%
% Calculate/keep quantile values in adaptive denoising.
%=========================================================================%

function qval = qf(df1, df2, alpha)
% Look-up table of F quantiles qf(1-alpha, df1, df2).
% df1  : 1..4   (numerator degrees of freedom)
% df2  : 1..40  (denominator degrees of freedom)
% alpha: significance level — one of 0.5, 0.1, 0.05, 0.025, 0.01
alpha_supported = [0.5, 0.1, 0.05, 0.025, 0.01];
aidx = find(abs(alpha_supported - alpha) < 1e-10);
if isempty(aidx)
    error('hmri_calc_paws:qf:unsupportedAlpha', ...
        'alpha must be one of: 0.5, 0.1, 0.05, 0.025, 0.01');
end

% quantile(df1, df2, alpha_idx)  —  values from R: qf(1-alpha, df1, 1:40)
quantile = zeros(4, 40, 5);

% alpha = 0.5  (50th percentile)
quantile(1,:,1) = [  1.000000   0.666667   0.585060   0.548632   0.528074   0.514890   0.505723   0.498982   0.493818   0.489737   0.486430   0.483696   0.481399   0.479441   0.477753   0.476283   0.474990   0.473845   0.472823   0.471906   0.471079   0.470329   0.469645   0.469019   0.468444   0.467915   0.467425   0.466971   0.466549   0.466155   0.465788   0.465443   0.465120   0.464816   0.464530   0.464259   0.464004   0.463762   0.463533   0.463315];
quantile(2,:,1) = [  1.500000   1.000000   0.881102   0.828427   0.798770   0.779763   0.766548   0.756828   0.749381   0.743492   0.738719   0.734772   0.731455   0.728627   0.726187   0.724062   0.722193   0.720538   0.719061   0.717735   0.716538   0.715452   0.714463   0.713557   0.712726   0.711959   0.711250   0.710593   0.709982   0.709412   0.708879   0.708381   0.707912   0.707472   0.707057   0.706666   0.706296   0.705946   0.705614   0.705298];
quantile(3,:,1) = [  1.709226   1.134943   1.000000   0.940534   0.907146   0.885783   0.870944   0.860038   0.851684   0.845081   0.839730   0.835306   0.831587   0.828418   0.825684   0.823303   0.821209   0.819354   0.817698   0.816213   0.814872   0.813655   0.812546   0.811532   0.810600   0.809741   0.808947   0.808210   0.807525   0.806887   0.806290   0.805731   0.805207   0.804713   0.804249   0.803810   0.803395   0.803003   0.802631   0.802278];
quantile(4,:,1) = [  1.822715   1.207107   1.063226   1.000000   0.964562   0.941913   0.926193   0.914645   0.905804   0.898817   0.893157   0.888478   0.884546   0.881195   0.878305   0.875787   0.873574   0.871612   0.869863   0.868293   0.866875   0.865589   0.864418   0.863345   0.862361   0.861453   0.860614   0.859835   0.859112   0.858437   0.857806   0.857216   0.856661   0.856140   0.855649   0.855185   0.854747   0.854332   0.853939   0.853566];

% alpha = 0.1  (90th percentile)
quantile(1,:,2) = [ 39.863458   8.526316   5.538319   4.544771   4.060420   3.775950   3.589428   3.457919   3.360303   3.285015   3.225202   3.176549   3.136205   3.102213   3.073186   3.048110   3.026232   3.006977   2.989900   2.974653   2.960956   2.948585   2.937356   2.927117   2.917745   2.909132   2.901192   2.893846   2.887033   2.880695   2.874784   2.869259   2.864083   2.859225   2.854655   2.850349   2.846285   2.842442   2.838804   2.835354];
quantile(2,:,2) = [ 49.500000   9.000000   5.462383   4.324555   3.779716   3.463304   3.257442   3.113118   3.006452   2.924466   2.859511   2.806796   2.763167   2.726468   2.695173   2.668171   2.644638   2.623947   2.605612   2.589254   2.574569   2.561314   2.549290   2.538332   2.528305   2.519096   2.510609   2.502761   2.495483   2.488716   2.482407   2.476512   2.470990   2.465809   2.460936   2.456346   2.452014   2.447920   2.444044   2.440369];
quantile(3,:,2) = [ 53.593245   9.161790   5.390773   4.190860   3.619477   3.288762   3.074072   2.923796   2.812863   2.727673   2.660229   2.605525   2.560273   2.522224   2.489788   2.461811   2.437434   2.416005   2.397022   2.380087   2.364888   2.351170   2.338727   2.327390   2.317017   2.307491   2.298712   2.290595   2.283069   2.276071   2.269548   2.263453   2.257744   2.252387   2.247350   2.242605   2.238128   2.233896   2.229890   2.226092];
quantile(4,:,2) = [ 55.832961   9.243416   5.342644   4.107250   3.520196   3.180763   2.960534   2.806426   2.692680   2.605336   2.536188   2.480102   2.433705   2.394692   2.361433   2.332745   2.307747   2.285772   2.266303   2.248934   2.233345   2.219274   2.206512   2.194882   2.184242   2.174469   2.165463   2.157136   2.149415   2.142235   2.135542   2.129288   2.123430   2.117934   2.112765   2.107896   2.103302   2.098959   2.094848   2.090950];

% alpha = 0.05  (95th percentile)
quantile(1,:,3) = [161.447639  18.512821  10.127964   7.708647   6.607891   5.987378   5.591448   5.317655   5.117355   4.964603   4.844336   4.747225   4.667193   4.600110   4.543077   4.493998   4.451322   4.413873   4.380750   4.351244   4.324794   4.300950   4.279344   4.259677   4.241699   4.225201   4.210008   4.195972   4.182964   4.170877   4.159615   4.149097   4.139252   4.130018   4.121338   4.113165   4.105456   4.098172   4.091279   4.084746];
quantile(2,:,3) = [199.500000  19.000000   9.552094   6.944272   5.786135   5.143253   4.737414   4.458970   4.256495   4.102821   3.982298   3.885294   3.805565   3.738892   3.682320   3.633723   3.591531   3.554557   3.521893   3.492828   3.466800   3.443357   3.422132   3.402826   3.385190   3.369016   3.354131   3.340386   3.327654   3.315830   3.304817   3.294537   3.284918   3.275898   3.267424   3.259446   3.251924   3.244818   3.238096   3.231727];
quantile(3,:,3) = [215.707345  19.164292   9.276628   6.591382   5.409451   4.757063   4.346831   4.066181   3.862548   3.708265   3.587434   3.490295   3.410534   3.343889   3.287382   3.238872   3.196777   3.159908   3.127350   3.098391   3.072467   3.049125   3.027998   3.008787   2.991241   2.975154   2.960351   2.946685   2.934030   2.922277   2.911334   2.901120   2.891564   2.882604   2.874187   2.866266   2.858796   2.851741   2.845068   2.838745];
quantile(4,:,3) = [224.583241  19.246794   9.117182   6.388233   5.192168   4.533677   4.120312   3.837853   3.633089   3.478050   3.356690   3.259167   3.179117   3.112250   3.055568   3.006917   2.964708   2.927744   2.895107   2.866081   2.840100   2.816708   2.795539   2.776289   2.758710   2.742594   2.727765   2.714076   2.701399   2.689628   2.678667   2.668437   2.658867   2.649894   2.641465   2.633532   2.626052   2.618988   2.612306   2.605975];

% alpha = 0.025  (97.5th percentile)
quantile(1,:,4) = [647.789011  38.506329  17.443443  12.217863  10.006982   8.813101   8.072669   7.570882   7.209283   6.936728   6.724130   6.553769   6.414254   6.297939   6.199501   6.115127   6.042013   5.978052   5.921631   5.871494   5.826648   5.786299   5.749805   5.716639   5.686366   5.658624   5.633109   5.609564   5.587768   5.567535   5.548702   5.531129   5.514693   5.499288   5.484820   5.471206   5.458372   5.446254   5.434793   5.423937];
quantile(2,:,4) = [799.500000  39.000000  16.044106  10.649111   8.433621   7.259856   6.541520   6.059467   5.714705   5.456396   5.255889   5.095867   4.965266   4.856698   4.765048   4.686665   4.618874   4.559672   4.507528   4.461255   4.419918   4.382768   4.349202   4.318726   4.290932   4.265483   4.242094   4.220525   4.200572   4.182061   4.164840   4.148779   4.133765   4.119700   4.106496   4.094076   4.082372   4.071326   4.060881   4.050992];
quantile(3,:,4) = [864.162972  39.165495  15.439182   9.979199   7.763589   6.598799   5.889819   5.415962   5.078119   4.825621   4.630025   4.474185   4.347178   4.241728   4.152804   4.076823   4.011163   3.953863   3.903428   3.858699   3.818761   3.782886   3.750486   3.721080   3.694273   3.669736   3.647192   3.626408   3.607187   3.589359   3.572778   3.557318   3.542868   3.529334   3.516631   3.504685   3.493429   3.482807   3.472766   3.463260];
quantile(4,:,4) = [899.583310  39.248418  15.100979   9.604530   7.387886   6.227161   5.522594   5.052632   4.718078   4.468342   4.275072   4.121209   3.995898   3.891914   3.804271   3.729417   3.664754   3.608344   3.558706   3.514695   3.475408   3.440126   3.408268   3.379359   3.353009   3.328894   3.306741   3.286321   3.267438   3.249925   3.233640   3.218456   3.204267   3.190977   3.178505   3.166777   3.155728   3.145301   3.135445   3.126114];

% alpha = 0.01  (99th percentile)
quantile(1,:,5) = [4052.180695  98.502513  34.116222  21.197690  16.258177  13.745023  12.246383  11.258624  10.561431  10.044289   9.646034   9.330212   9.073806   8.861593   8.683117   8.530965   8.399740   8.285420   8.184947   8.095958   8.016597   7.945386   7.881134   7.822871   7.769798   7.721254   7.676684   7.635619   7.597663   7.562476   7.529766   7.499281   7.470801   7.444136   7.419117   7.395597   7.373445   7.352545   7.332794   7.314100];
quantile(2,:,5) = [4999.500000  99.000000  30.816520  18.000000  13.273934  10.924767   9.546578   8.649111   8.021517   7.559432   7.205713   6.926608   6.700965   6.514884   6.358873   6.226235   6.112114   6.012905   5.925879   5.848932   5.780416   5.719022   5.663699   5.613591   5.567997   5.526335   5.488118   5.452937   5.420445   5.390346   5.362385   5.336343   5.312029   5.289277   5.267941   5.247894   5.229022   5.211225   5.194413   5.178508];
quantile(3,:,5) = [5403.352014  99.166201  29.456695  16.694369  12.059954   9.779538   8.451285   7.590992   6.991917   6.552313   6.216730   5.952545   5.739380   5.563886   5.416965   5.292214   5.185000   5.091890   5.010287   4.938193   4.874046   4.816606   4.764877   4.718051   4.675465   4.636570   4.600907   4.568091   4.537795   4.509740   4.483686   4.459429   4.436787   4.415606   4.395749   4.377096   4.359540   4.342988   4.327356   4.312569];
quantile(4,:,5) = [5624.583330  99.249372  28.709898  15.977025  11.391928   9.148301   7.846645   7.006077   6.422085   5.994339   5.668300   5.411951   5.205330   5.035378   4.893210   4.772578   4.668968   4.579036   4.500258   4.430690   4.368815   4.313429   4.263567   4.218445   4.177420   4.139960   4.105622   4.074032   4.044873   4.017877   3.992811   3.969477   3.947701   3.927333   3.908241   3.890308   3.873433   3.857524   3.842502   3.828294];

qval = quantile(df1, df2, aidx);
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

%% =======================================================================%
% modify get_trtefa to get TE values
% To extract the TE values out of the descriptor field of 
% the nifti header or the extended header/JSON metadata.
%=========================================================================%
function p = get_te(P)

N = nifti(P);
nN = numel(N);
p(nN) = struct('te',[]);

for ii = 1:numel(N)
    p(ii).te = get_metadata_val(P(ii,:),'EchoTime');
end

end
