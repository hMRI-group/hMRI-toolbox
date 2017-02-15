function isHdrExtended = has_extended_header(niifilenam)
% Must strip the ',1' (at the end of the file extension '.nii,1') 
% if the name has been collected using e.g. spm_select:
niifilenam = spm_file(niifilenam,'number','');
% Check whether a NIFTI file has extended header or not...
fid = fopen(niifilenam,'r+');
std_hdr_size = 348;
fseek(fid, std_hdr_size, 'bof');
isHdrExtended = fread(fid,1,'uint8');
fclose(fid);