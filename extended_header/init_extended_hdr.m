function N = init_extended_hdr(N, hdr)
% To create the JSON-encoded extended header during DICOM to nifti
% conversion, including all acquisition parameters. This function is
% called by spm_dicom_convert and is applied to the nifti object being
% created. The DICOM header returned by spm_dicom_headers is tidied up and
% converted into JSON. The size of the JSON header is used to define a new
% offset (N.dat.offset) for the image data in the nifti file and the JSON
% header is written into the nifti file. The modified nifti object N is
% returned so the data can be written in it according to the new offset (in
% spm_dicom_convert).
%__________________________________________________________________________
% FORMAT N = init_extended_hdr(N, hdr)
% hdr       a single matlab structure containing the header (from
%           spm_dicom_headers)
% N(input)  the nifti object created in spm_dicom_convert with file name
%           and default offset
% N(input)  the nifti object modified with extended header and
%           corresponding extended offset
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

% check JSONLAB is present
if ~exist('savejson.m','file')
    warning('JSONLAB does not seem to be in the Matlab Path.');
    fprintf('JSONLAB is required to set and get extended headers.\nThe extended header will not be created.\n');
    %jhdr = '';
    %offset = 352;
    return;
end

% TIDY UP THE DICOM HEADER:
% 1) CSA fields
spm_tdyhdr = eb_spm_tidycsa(hdr);
% 2) MrPhoenixProtocol field and ASCCONV part of the DICOM header
if isfield(spm_tdyhdr,'CSASeriesHeaderInfo')
    if isfield(spm_tdyhdr.CSASeriesHeaderInfo,'MrPhoenixProtocol')
        spm_tdyhdr.CSASeriesHeaderInfo.MrPhoenixProtocol = eb_read_phoenix(spm_tdyhdr.CSASeriesHeaderInfo.MrPhoenixProtocol);
    end
end
% 3) basic anonymisation
if isfield(spm_tdyhdr,'PatientName')
    spm_tdyhdr.PatientName = 'anonymous';
end
if isfield(spm_tdyhdr,'PatientBirthDate')
    %t1 = datenum(spm_tdyhdr.PatientBirthDate,'yyyymmdd');
    %t2 = datenum(spm_tdyhdr.StudyDate,'yyyymmdd');
    %spm_tdyhdr.PatientAge = round((t2-t1)*10/365.25)/10;
    spm_tdyhdr = rmfield(spm_tdyhdr,'PatientBirthDate');
end

% ORGANIZE STRUCTURE WITH acqpar AND history FIELDS
exthdr.history.step{1} = struct('descrip','dicom to nifti import', ...
                                'version','v0.0', ...
                                'procpar',[], ...
                                'imtype','unprocessed image', ...
                                'units','a.u.');
exthdr.acqpar = spm_tdyhdr;

% JSONify the header
jhdr = savejson('',exthdr);

% only bother with the following if the jhdr isn't empty
if ~isempty(jhdr);
    % a few parameters:
    % standard header is always 348 bytes long:
    std_hdr_size = 348;
    % the 'extension' array field (bytes 348-351)
    isHdrExtended = [1 0 0 0]; % [0 0 0 0] if no extended header is present
    % the extension includes the json header + two 32-bit numbers = 8 bytes
    % (the size of the extension and the ID of the extension):
    ext_hdr_size = length(jhdr)+8;
    % ID of the extension (32-bit number) = anything > 4, arbitrarily set to 27
    % for now, just because I like this number :)...
    ext_hdr_id = 27;
    % the offset must be >352+ext_hdr_size and a multiple of 16
    offset = ceil((352+ext_hdr_size)/16)*16;
    
    % modify the offset to write data in the nifti file
    N.dat.offset = offset;
    % since offset has been modified, N must be re-created:
    create(N);
    
    % open nifti file to write the extended header
    fid = fopen(N.dat.fname,'r+');
    % standard header is always 348 bytes long, move first there:
    fseek(fid, std_hdr_size, 'bof');
    % the next 4 bytes indicate whether there is an extension.
    fwrite(fid,isHdrExtended,'uint8');
    % we should now be @byte 348+4 = 352: ftell(fid)
    % write the 32-bit numbers
    fwrite(fid,ext_hdr_size,'uint32');
    fwrite(fid,ext_hdr_id,'uint32');
    % we should now be @byte 348+4+8 = 360: ftell(fid)
    % write the jhdr
    fprintf(fid,'%s',jhdr); % disp(jhdr);
    fclose(fid);
end

