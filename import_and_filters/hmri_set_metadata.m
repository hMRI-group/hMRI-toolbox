function hmri_set_metadata(filelist, mstruc, varargin)
% To insert or update JSON-encoded metadata into nifti images (for extended
% nifti format) or JSON files associated to standard nifti format.
%__________________________________________________________________________
% FORMAT hmri_set_metadata(filelist, mstruc[, mtype])
% filelist  the name of a nifti image file, or list of files
% mstruc    a single matlab structure containing the metadata (new or
%           updated). The same metadata will be fed into each file.
% mtype     type of metadata storage used (either 'json' for separate json
%           file or 'nii+' for extended nifti header). Default value is 
%           'json' if the input filename has json extension or 'nii+' if
%           the input filename has nii extension.
%__________________________________________________________________________
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
% - Enrico Reimer, Max Planck Institute for Human Cognitive and Brain
%   Sciences, Leipzig, Germany 
%__________________________________________________________________________

% Must strip the ',1' (at the end of the file extension '.nii,1') 
% if the name has been collected using e.g. spm_select:
filelist = spm_file(filelist,'number','');

% Set compacting options for hmri_jsonwrite (not compact for better
% readability):
opts.compact = false;

for cfile = 1:size(filelist,1)
    [pth, fnam, ext] = fileparts(filelist(cfile,:));
    
    % define mtype for the current file and check whether mtype and input
    % file extension are compatible...
    mtype = 'json';
    if nargin==2
        switch lower(ext)
            case '.json'
                mtype = 'json';
            case '.nii'
                mtype = 'nii+';
        end
    else
        mtype = varargin{1};
        if strcmp(mtype,'nii+') && strcmp(ext,'.json')
            error('Metadata cannot be added as extended header to a JSON file!');
        end
    end
    
%     % check whether there is already existing metadata available and
%     % whether they can be overwritten
%     if exist(filelist(cfile,:),'file')
%         oldmstruc = hmri_get_metadata(filelist(cfile,:));
%         if ~isempty(oldmstruc{1})
%             warning('Existing metadata already stored in %s or associated JSON file. This operation will overwrite the existing metadata.', filelist(cfile,:));
%             ok = input('Are you sure you want to proceed? (Y/N) ','s');
%             if lower(ok)=='n', return; end
%         end
%     end
    
    switch mtype
        case 'json'
            switch lower(ext)
                case '.nii'
                    hmri_jsonwrite(fullfile(pth,[fnam '.json']),mstruc,opts);
                case '.json'
                    hmri_jsonwrite(filelist(cfile,:),mstruc,opts);
                otherwise
                    error('Unknown file extension (%s). JSON metadata can be stored as nifti extended header (*.nii) or JSON file (*.json)',ext);
            end
            
        case 'nii+'
            % JSONify the matlab metadata structure
            jsonhdr = hmri_jsonwrite(mstruc,opts);

            % set a few parameters:
            % standard header is always 348 bytes long:
            std_hdr_size = 348;
            % the 'extension' array field (bytes 348-351)
            isHdrExtended = [1 0 0 0];
            % the length of the extension includes two 32-bit numbers = 8 bytes (the
            % size itself and the ID of the extension):
            ext_hdr_size = length(jsonhdr)+8;
            % ID of the extension (32-bit number) = anything > 4, arbitrarily set to 27
            % for now, just because I like this number :)...
            ext_hdr_id = 27;

            % read and save the data for later
            Y = spm_read_vols(spm_vol(filelist(cfile,:)));
            % modify the nifti file offset according to the length of the extension
            N = nifti(filelist(cfile,:));
            % the offset must be >352+ext_hdr_size and a multiple of 16
            N.dat.offset = ceil((352+ext_hdr_size)/16)*16;
            % write the updated nifti standard header
            create(N);
            % rewrite the data (according to the updated offset)
            N.dat(:,:,:) = Y;
            % open nifti file to write the extended header
            fid = fopen(filelist(cfile,:),'r+');
            % standard header is always 348 bytes long, move first there:
            fseek(fid, std_hdr_size, 'bof');
            % the next 4 bytes indicate whether there is an extension.
            fwrite(fid,isHdrExtended,'uint8');
            % we should now be @byte 348+4 = 352: ftell(fid)
            % write the 32-bit numbers
            fwrite(fid,ext_hdr_size,'uint32');
            fwrite(fid,ext_hdr_id,'uint32');
            % we should now be @byte 348+4+8 = 360: ftell(fid)
            % write the jsonhdr
            fprintf(fid,'%s',jsonhdr); % disp(jsonhdr);
            fclose(fid);
    end
end

end