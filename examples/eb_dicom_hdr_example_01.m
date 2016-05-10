function eb_dicom_hdr_example_01

% Trying to sort out a header extractor for DICOM conversion...
% This first example gather almost all the available information from the
% header and store it (after a minimum of trimming and tidying up) as a
% structure. The structure is then JSONified (since the metadata extended
% header will use JSON) and saved in a txt file. The length of the JSON
% string (in bytes) is printed within the txt file and should correspond to
% the txt file size. For convenience in this sandbox/testing/try&fail
% stage, the structure is also saved as a Matlab mat file. 

addpath(fullfile(pwd, '..', 'read_dicom_header'));

filename{1} = fullfile(pwd,'FDL_40_WILLIAM_THACKERAY_MR.MR.CRC_PROTOCOLS.0002.0001.2016.02.03.13.39.48.883561.36253775.IMA');
filename{2} = fullfile(pwd,'FDL_40_WILLIAM_THACKERAY_MR.MR.CRC_PROTOCOLS.0008.0091.2016.02.03.13.39.48.883561.36410704.IMA');
filename{3} = fullfile(pwd,'FDL_40_WILLIAM_THACKERAY_MR.MR.CRC_PROTOCOLS.0009.0020.2016.02.03.13.39.48.883561.36430073.IMA');
filename{4} = fullfile(pwd,'DUMANOIR_PHANTOM.MR.CRC_PHYSICS.0002.0004.2015.07.01.10.57.37.800923.67595629.IMA');

for i=1%:length(filename)
    % read header with SPM tool
    spm_hdr = spm_dicom_headers(filename{i});
    % tidy it up:
    % 1) CSA fields
    spm_hdr = eb_spm_tidycsa(spm_hdr{1});
    % 2) MrPhoenixProtocol field and ASCCONV part of the DICOM header
    spm_hdr.CSASeriesHeaderInfo.MrPhoenixProtocol = eb_read_phoenix(spm_hdr.CSASeriesHeaderInfo.MrPhoenixProtocol);
    
    % JSONify
    json_spm_hdr = savejson('',spm_hdr);
    [pth,nam,~] = fileparts(filename{i});
        
    fid = fopen(fullfile(pth,[nam '_spmtdyhdr.txt']),'w');
    fprintf(fid,'JSON header from spm_dicom_headers (length = %d)\n',length(json_spm_hdr));
    fprintf(fid,'%s',json_spm_hdr);
    fclose(fid);
    
    save(fullfile(pth,[nam '_spmtdyhdr.mat']), 'spm_hdr');
end