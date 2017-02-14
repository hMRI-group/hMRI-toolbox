function set_metadata(filelist, mstruc, varargin)
% To insert or update JSON-encoded metadata into nifti images (for extended
% nifti format) or JSON files associated to standard nifti format. 
%==========================================================================
% FORMAT hmri_set_metadata(filelist, mstruc[, json])
% filelist  the name of a nifti image file, a json metadata file, or list
%           of such files (can be mixed, but...). 
% mstruc    a single matlab structure containing the metadata (new or
%           updated). The same metadata will be fed into each file.
% json      a structure with fields
%               extended: true/false -> JSON metadata are stored as
%                           extended header in the nii file. Default is
%                           true if input file has nii extension, false
%                           otherwise. 
%               separate: true/false -> JSON metadata are stored as a
%                           separate JSON file. Default is true if input
%                           file has json extension, false otherwise. 
%               overwrite: true/false -> force overwrite existing metadata
%                           if true, keep existing metadata and add new
%                           ones if false. Default is false. However, it 
%                           is recommended to use json.overwrite = true 
%                           and update the metadata in a more readable way
%                           (good housekeeping recommendation - see manual
%                           about metadata structure!). The proposed
%                           non-overwritng default is intended to prevent
%                           metadata loss if overwritten by mistake... 
%
%==========================================================================
% NOTES: 
% If no json options are provided, the script will search for nii file and
% associated json file and update both the extended nii header (if any) and
% the json associated file (if any). If no metadata yet, by default, a json
% file is created (with the file name provided and json extension).
%
% If json options are provided, the script will update/create the metadata
% as specified and throw a warning if existing metadata, not specified in
% the json options, are present but not specified (hence not updated). For
% example, if a json file exists but json.separate = false, the json
% metadata won't be updated and a warning will be issued.  
%
%==========================================================================
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
% - Enrico Reimer, Max Planck Institute for Human Cognitive and Brain
%   Sciences, Leipzig, Germany 
%==========================================================================

% Must strip the ',1' (at the end of the file extension '.nii,1') 
% if the name has been collected using e.g. spm_select:
filelist = spm_file(filelist,'number','');

% Set indent options for spm_jsonwrite (using tabs and non-compact
% formatting for better readability): 
opts.indent = '\t';

for cfile = 1:size(filelist,1)
    [pth, fnam, ~] = fileparts(filelist(cfile,:));
    
    % check status for current file (exisiting json/nii files, nii file
    % with extended or not extended header...)
    jsonexist = 0;niiexist = 0; niiextended = 0;
    if exist(fullfile(pth, [fnam '.json']),'file'); jsonexist = 1;end
    if exist(fullfile(pth, [fnam '.nii']),'file')
        niiexist = 1;
        % check whether extended nii format
        niiextended = has_extended_header(fullfile(pth, [fnam '.nii']));
    end
    
    % define json structure for the current file and check for conflicts
    % with existing metadata... 
    json = struct('extended',false,'separate',false,'overwrite',false);
    if nargin==2
        if jsonexist; json.separate = true; end
        if niiextended; json.extended = true; end
        if ~jsonexist && ~niiextended; json.separate = true; end
    else
        if isfield(varargin{1},'extended'); json.extended = varargin{1}.extended; end
        if isfield(varargin{1},'separate'); json.separate = varargin{1}.separate; end
        if isfield(varargin{1},'overwrite'); json.overwrite = varargin{1}.overwrite; end
        if ~json.extended && niiextended
            warning('Existing metadata in the nifti extended header won''t be updated');
        end
        if ~json.separate && jsonexist
            warning('Existing metadata in the associated JSON file won''t be updated');
        end
        if ~json.separate && ~json.extended
            warning('JSON options not properly set. Metadata will be saved/updated in a separate JSON file.');
            json.separate = true;
        end
        if json.extended && ~niiexist
            warning('There is no nii file available to save metadata as extended NIFTI header. Metadata will be saved/updated in a separate JSON file');
            json.separate = true;
            json.extended = false;
        end
    end
    
    % deal with overwriting specifications...
    if ~json.overwrite
        existing_mstruc = get_metadata(filelist(cfile,:));
        if ~isempty(existing_mstruc{1})
            n=1; archivefieldname = sprintf('PreviousMetadata%.3d',n);
            while isfield(mstruc,archivefieldname)
                n=n+1; archivefieldname = sprintf('PreviousMetadata%.3d',n);
            end
            mstruc.(archivefieldname) = existing_mstruc{1};
        end
    end
    
    % pfew, now let's do the job...
    if json.separate
        spm_jsonwrite(fullfile(pth, [fnam '.json']),mstruc,opts);
    end
    if json.extended
        % read and store the data for later
        Y = spm_read_vols(spm_vol(fullfile(pth, [fnam '.nii'])));
        % create handle to NIFTI object
        N = nifti(fullfile(pth, [fnam '.nii']));
        % create JSONified header and calculate NIFTI extended offset
        [jhdr, offset] = get_jhdr_and_offset(mstruc);
        % modify the nifti file offset 
        N.dat.offset = offset;
        % make this change effective by rewriting the NIFTI header
        create(N);
        % rewrite the data (done according to the updated offset)
        for cv = 1:size(Y,4)
            N.dat(:,:,:,cv) = Y(:,:,:,cv);
        end
        % write the extended header
        write_extended_header(fullfile(pth, [fnam '.nii']),jhdr);
    end
end