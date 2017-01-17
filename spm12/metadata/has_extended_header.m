function isHdrExtended = has_extended_header(niifilenam)
% Check whether a NIFTI file has extended header or not...
fid = fopen(niifilenam,'r+');
std_hdr_size = 348;
fseek(fid, std_hdr_size, 'bof');
isHdrExtended = fread(fid,1,'uint8');
fclose(fid);