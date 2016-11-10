function hMRI_set_extended_hdr(filelist, hdr)
% To insert or update JSON-encoded metadata into nifti images.
%__________________________________________________________________________
% FORMAT set_extended_hdr(filelist, hdr)
% filelist    the name of a nifti image file, or list of files
% hdr         a single matlab structure containing the metadata (new or
%               updated). The same metadata will be fed into each file.
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

% JSONify the header
jsonhdr = savejson('',hdr);

% few parameters:
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

for cfile = 1:size(filelist,1)
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