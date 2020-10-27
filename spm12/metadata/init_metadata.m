function N = init_metadata(N, hdr, json)
% To create JSON-encoded metadata during DICOM to nifti conversion,
% including all acquisition parameters. The metadata can either (or both)
% be stored in the extended header of the nifti image or saved as a
% separate file. This function is called by spm_dicom_convert. In case of
% an extended nii header, the size of the JSON header is used to define a
% new offset (N.dat.offset) for the image data in the nifti file and the
% JSON header is written into the nifti file. The modified nifti object N
% is returned so the data can be written in it according to the new offset
% (in spm_dicom_convert). In case of a separate JSON file, N is returned
% unchanged and the JSON metadata are written in a separate file (same file
% name as the nifti image, with .json extension).
%__________________________________________________________________________
% FORMAT N = init_metadata(N, hdr, json)
% hdr       a single matlab structure containing the header (from
%           spm_dicom_headers)
% N(input)  the nifti object created in spm_dicom_convert with file name
%           and default offset
% json      a structure with fields
%               extended: true/false -> JSON metadata are stored as
%                           extended header in the nii file.
%               separate: true/false -> JSON metadata are stored as a
%                           separate JSON file.
%               anonym: 'basic','full','none' (see anonymise_metadata.m)
% N(output) the nifti object modified with extended header and
%           corresponding extended offset (if json.extended = true).
%           N(output) = N(input) if json.extended = false.
%==========================================================================
% Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
%==========================================================================

% INITIALIZE AND ORGANIZE STRUCTURE WITH acqpar AND history FIELDS
dicom_convert_version = {sprintf('%s %s', spm_check_version, version),sprintf('spm_dicom_convert.m - %s', spm('Version'))};
metadata.history.procstep = struct('descrip','dicom to nifti import', 'version', {dicom_convert_version}, 'procpar', []);
metadata.history.input(1) = struct('filename','AnonymousFileName', 'history',[]);
if isfield(hdr,'ImageType')
    metadata.history.output = struct('imtype',hdr.ImageType, 'units','a.u.');
else
    metadata.history.output = struct('imtype','Unprocessed MR image', 'units','a.u.');
end

% Tidy up and rearrange CSA fields in the header, including formatting the
% ASCII part into proper Matlab structure (Note: this is specific to
% Siemens DICOM format)
hdr = reformat_spm_dicom_header(hdr);

% NB: Anonymisation is only basic and might not be effective!!!
hdr = anonymise_metadata(hdr,struct('anonym',json.anonym));
metadata.acqpar = hdr;

if json.extended
    % JSONify the header and calculate the required offset
    [jhdr, jhdr_size, offset] = get_jhdr_and_offset(metadata);
    % modify the offset of the nifti
    N.dat.offset = offset;
    % make the offset modification effective by rewriting the standard
    % nifti header, including the offset:
    create(N);
    % write the extended header (can be done before data are written or
    % after, does not matter:
    write_extended_header(N.dat.fname,jhdr,jhdr_size);
end

if json.separate
    [pth,fnam,ext] = fileparts(N.dat.fname);
    spm_jsonwrite(fullfile(pth,[fnam '.json']),metadata, struct('indent','\t'));
end
