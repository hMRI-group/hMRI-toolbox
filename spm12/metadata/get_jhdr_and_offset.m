function [jhdr, offset] = get_jhdr_and_offset(hdr)
%==========================================================================
% FORMAT
% [jhdr, offset] = get_jhdr_and_offset(hdr)
% hdr       a Matlab structure
% jhdr      the JSONified Matlab structure
% offset    the offset at which image data start in the nifti file
%==========================================================================
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
%==========================================================================

% JSONify the Matlab structure
jhdr = spm_jsonwrite(hdr,struct('indent','\t'));

% the length of the extension includes two 32-bit numbers = 8 bytes (the
% size of the extension itself and the ID of the extension):
ext_hdr_size = length(jhdr)+8;

% the offset must be >352+ext_hdr_size AND a multiple of 16
offset = ceil((352+ext_hdr_size)/16)*16;

% we fill up the jsonhdr with white spaces to avoid confusion when
% human-reading the extended header
needed_white_spaces = offset - (352+ext_hdr_size);
jhdr = [jhdr repmat(' ',1,needed_white_spaces)];