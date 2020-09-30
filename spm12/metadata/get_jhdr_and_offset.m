function [jhdr, jhdr_size, offset] = get_jhdr_and_offset(hdr)
%==========================================================================
% FORMAT
% [jhdr, jhdr_size, offset] = get_jhdr_and_offset(hdr)
% hdr       a Matlab structure
% jhdr      the JSONified Matlab structure
% jhdr_size size of the JSON header encoding in bytes
% offset    the offset at which image data start in the nifti file
%==========================================================================
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
%==========================================================================

% JSONify the Matlab structure
jhdr = spm_jsonwrite(hdr,struct('indent','\t'));

% simulate header write, fprintf() returns actual bytes written. length()
% doesn't work here because characters might be larger than one byte (UTF-8)
sim_fname = tempname();
fid = fopen(sim_fname, 'w');
jhdr_size = fprintf(fid, '%s', jhdr);
fclose(fid);
delete(sim_fname);

% the length of the extension includes two 32-bit numbers = 8 bytes (the
% size of the extension itself and the ID of the extension):
ext_hdr_size = jhdr_size+8;

% the offset must be >352+ext_hdr_size AND a multiple of 16
offset = ceil((352+ext_hdr_size)/16)*16;

% we fill up the jsonhdr with white spaces to avoid confusion when
% human-reading the extended header
needed_white_spaces = offset - (352+ext_hdr_size);
jhdr = [jhdr repmat(' ',1,needed_white_spaces)];
