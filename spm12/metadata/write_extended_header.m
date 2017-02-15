function write_extended_header(fnam, jhdr)
%==========================================================================
% FORMAT
% write_extended_header(fnam, jsonhdr)
% fnam      the name of a nifti file
% jsonhdr   the JSONified Matlab structure to be written as extended header
%==========================================================================
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
%==========================================================================

% set a few parameters for NIFTI format:
% standard header is always 348 bytes long:
std_hdr_size = 348;
% the 'extension' array field (bytes 348-351)
isHdrExtended = [1 0 0 0];
% the length of the extension includes two 32-bit numbers = 8 bytes (the
% size itself and the ID of the extension):
ext_hdr_size = length(jhdr)+8;
% ID of the extension (32-bit number) = anything > 4, arbitrarily set to 27
% for now, just because I like this number :)...
ext_hdr_id = 27;

% open nifti file to write the extended header
fid = fopen(fnam,'r+');
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
fprintf(fid,'%s',jhdr); % disp(jsonhdr);
fclose(fid);