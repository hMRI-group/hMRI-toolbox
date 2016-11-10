function hdr = hMRI_get_extended_hdr(filelist)
% To retrieve JSON-encoded metadata from nifti images.
%__________________________________________________________________________
% FORMAT hdr = get_extended_hdr(filelist)
% filelist    the name of a nifti image file, or list of files
% hdr         a cell array of matlab structures
%__________________________________________________________________________
% DEPENDENCIES
% - JSONlab
%       An open-source MATLAB/Octave JSON/UBJSON encoder and decoder
%       http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?jsonlab
%       https://github.com/fangq/jsonlab.git
%       Copyright (C) 2011-2015 Qianqian Fang <fangq@nmr.mgh.harvard.edu>
%       License: BSD or GNU General Public License version 3 (GPL v3),
%       see License*.txt 
%       Version: 1.2 (Optimus - Update 2)
%__________________________________________________________________________
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
% - Enrico Reimer, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig, Germany
%__________________________________________________________________________

% Must strip the ',1' (at the end of the file extension '.nii,1') 
% if the name has been collected using e.g. spm_select:
filelist = spm_file(filelist,'number','');

for cfile = 1:size(filelist,1)
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
        fprintf('\nNifti header extension available in file \n%s \nisHdrExtended = %d\n',filelist(cfile,:),isHdrExtended(1));
        % status = fseek(fid, std_hdr_size+4, 'bof'); % not necessary since we've already moved at the right position while reading
        ext_hdr_size = fread(fid,1,'uint32'); % read size of the extension (32-bit number) = length of the JSON string
        ext_hdr_id = fread(fid,1,'uint32'); % read the ID of the extension (32-bit number) = anything > 4, arbitraryly set to 27 for now, just because I like this number :)...
        jsonhdr = fscanf(fid,'%c',ext_hdr_size-8); % disp(jsonhdr);
        hdr{cfile} = loadjson(jsonhdr); % disp(hdr);
    else
        fprintf('\nNo nifti header extension available in file \n%s \nisHdrExtended = %d\n',filelist(cfile,:),isHdrExtended(1));
        hdr{cfile} = [];
    end
    fclose(fid);
end

end