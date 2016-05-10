function [jhdr, offset] = init_extended_hdr(hdr)
% To create the JSON-encoded extended header during DICOM to nifti
% convertions, including all acquisition parameters. This function is
% called by spm_dicom_convert and returne the JSONified header and
% corresponding offset for insertion as extended header.
%__________________________________________________________________________
% FORMAT [jhdr, offset] = init_extended_hdr(hdr)
% hdr       a single matlab structure containing the header (from
%           spm_dicom_headers) 
% jhdr      the tidied up, JSONified header to be inserted as extended header 
% offset    the corresponding offset to write data in the nifti image
%__________________________________________________________________________
% DEPENDENCIES
% - JSONlab
%       An open-source MATLAB/Octave JSON/UBJSON encoder and decoder
%       http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?jsonlab
%       https://github.com/fangq/jsonlab.git
%       Copyright (C) 2011-2015 Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%       License: BSD or GNU General Public License version 3 (GPL v3), see License*.txt
%       Version: 1.2 (Optimus - Update 2)
% - scripts to read and tidy up the DICOM header (in read_dicom_header)
%       - eb_read_ASCII
%       - eb_read_phoenix
%       - eb_spm_tidycsa
%       - get_numaris4_val (from spm_dicom_convert)
%__________________________________________________________________________
% Written by
% - Evelyne Balteau, Cyclotron Research Centre, Liège, Belgium
% - Enrico Reimer, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig, Germany
%__________________________________________________________________________

% tidy up the header:
% 1) CSA fields
spm_tdyhdr = eb_spm_tidycsa(hdr);
% 2) MrPhoenixProtocol field and ASCCONV part of the DICOM header
if isfield(spm_tdyhdr,'CSASeriesHeaderInfo')
    if isfield(spm_tdyhdr.CSASeriesHeaderInfo,'MrPhoenixProtocol')
        spm_tdyhdr.CSASeriesHeaderInfo.MrPhoenixProtocol = eb_read_phoenix(spm_tdyhdr.CSASeriesHeaderInfo.MrPhoenixProtocol);
    end
end

% check JSONLAB is present
if ~exist('savejson')
    warning('JSONLAB does not seem to be in the Matlab Path.');
    fprintf('JSONLAB is required to set and get extended headers.\nThe extended header will not be created.\n');
    jhdr = '';
    offset = 352;
else
    % JSONify the header
    jhdr = savejson('',spm_tdyhdr);

    % few parameters:
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
    % the offset must be >352+ext_hdr_size and a multiple of 16
    offset = ceil((352+ext_hdr_size)/16)*16;
    
end