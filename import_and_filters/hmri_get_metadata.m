function hdr = hmri_get_metadata(varargin)
% PURPOSE
% To retrieve JSON-encoded metadata from extended nifti images or
% associated JSON file.
%
% FORMAT 
% hdr = hmri_get_metadata(filelist)
% filelist  the name of a nifti image file (extended nifti or simple nifti
%           accompanied by JSON file with identical file name), associated
%           JSON file, or list of such files (can be mixed list). If no
%           argument passed, files can be selected using spm_select.
% hdr       a cell array of matlab structures
%
% DEPENDENCIES
% spm_jsonread
%--------------------------------------------------------------------------
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
% - Enrico Reimer, Max Planck Institute for Human Cognitive and Brain
%   Sciences, Leipzig, Germany 

if nargin==1
    filelist = varargin{1};
elseif nargin==0
    filelist = spm_select(Inf,'any','Select files...');
end

% Must strip the ',1' (at the end of the file extension '.nii,1') 
% if the name has been collected using e.g. spm_select:
filelist = spm_file(filelist,'number','');
hdr = cell(1,size(filelist,1));

for cfile = 1:size(filelist,1)
    
    % different behaviours according to metadata format (extended nii,
    % normal nii (no metadata) or separate JSON file)
    [pth,fnam,ext] = fileparts(filelist(cfile,:));
    switch lower(ext)
        case '.nii'
            % open nifti file
            fid = fopen(filelist(cfile,:),'r+');
            % standard header is always 348 bytes long, move first there:
            std_hdr_size = 348;
            fseek(fid, std_hdr_size, 'bof');
            % the next 4 bytes indicate whether there is an extension
            % (actually, only the first byte is of interest)
            isHdrExtended = fread(fid,4,'uint8');
            % we should now be @byte 348+4 = 352
            % fprintf('\nCurrent position in the file should be 348+4 = 352: ftell(fid) = %d\n',ftell(fid));
            
            if isHdrExtended(1)
                % fprintf('\nNifti header extension available in file \n%s \nisHdrExtended = %d\n',filelist(cfile,:),isHdrExtended(1));
                fprintf('\nNifti header extension available in file \n%s \n',filelist(cfile,:));
                % status = fseek(fid, std_hdr_size+4, 'bof'); % not necessary since we've already moved at the right position while reading
                ext_hdr_size = fread(fid,1,'uint32'); % read size of the extension (32-bit number) = length of the JSON string
                ext_hdr_id = fread(fid,1,'uint32'); %#ok<NASGU> % read the ID of the extension (32-bit number) = anything > 4, arbitraryly set to 27 for now, just because I like this number :)...
                jsonhdr = fscanf(fid,'%c',ext_hdr_size-8); % disp(jsonhdr);
                % convert the jsonhdr (string containing the JSON
                % structure) into a matlab structure using spm_jsonread:
                hdr{cfile} = spm_jsonread(jsonhdr);
            else
                % check whether there is a JSON file associated to the
                % non-extended nii image:
                if exist(fullfile(pth,[fnam '.json']),'file')
                    tmp = hmri_get_metadata(fullfile(pth,[fnam '.json']));
                    hdr{cfile} = tmp{1};
                else
                    fprintf(['\nNo nifti header extension available in file \n%s' ...
                        '\nNo associated JSON metadata file either.\n'],filelist(cfile,:));
                    hdr{cfile} = [];
                end
            end
            fclose(fid);
            
        case '.json'
            % can simply use spm_jsonread to convert the JSON structure
            % contained in the file into a matlab structure:
            try
                hdr{cfile} = spm_jsonread(filelist(cfile,:));
            catch MExc
                fprintf(1,'\n%s\n', MExc.getReport);
                fprintf('\nError when reading the JSON file \n%s \nFile might be empty or corrupted.\n',filelist(cfile,:));
                hdr{cfile} = [];
            end
            
        otherwise
            fprintf(['\nUnknown file format for file \n%s\n'...
                'The input file must be either nifti (preferably with associated JSON file),\n'...
                'extended nifti or JSON.\n'], filelist(cfile,:));
            hdr{cfile} = [];
    end
end

end