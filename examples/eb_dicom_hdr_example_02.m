function eb_dicom_hdr_example_02

% Converting a simple EPI volume into nifti with (very much) extended
% metadata header.

addpath(fullfile(pwd, '..'));
addpath(fullfile(pwd, '..','read_dicom_header'));

filename{1} = fullfile(pwd,'DUMANOIR_PHANTOM.MR.CRC_PHYSICS.0002.0004.2015.07.01.10.57.37.800923.67595629.IMA');

for i=1:length(filename)
    % read header with SPM tool
    spm_hdr = spm_dicom_headers(filename{i});
    % tidy it up:
    % 1) CSA fields
    spm_tdyhdr = eb_spm_tidycsa(spm_hdr{1});
    % 2) MrPhoenixProtocol field and ASCCONV part of the DICOM header
    spm_tdyhdr.CSASeriesHeaderInfo.MrPhoenixProtocol = eb_read_phoenix(spm_tdyhdr.CSASeriesHeaderInfo.MrPhoenixProtocol);

    % now convert IMA file into nifti
    out = spm_dicom_convert(spm_hdr,'mosaic','flat','nii');

    % and insert extended header
    set_extended_hdr(out.files{1},spm_tdyhdr);
end