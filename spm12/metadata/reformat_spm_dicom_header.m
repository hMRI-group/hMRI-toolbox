function hdr = reformat_spm_dicom_header(hdr)
%-------------------------------------------------------------------------%
% To tidy up and rearrange CSA fields in the header, including formatting
% the ASCII part into a proper Matlab structure (Note: this is specific to
% Siemens DICOM format but could be extended to other cases where needed) 
%-------------------------------------------------------------------------%
% FORMAT hdrout = reformat_spm_dicom_header(hdrin)
% where hdr is a dicom header structure as output by spm_dicom_headers.
%-------------------------------------------------------------------------%

% list of CSA-like fields:
CSAlist = {'CSAImageHeaderInfo', 'CSASeriesHeaderInfo','CSANonImageHeaderInfoVA','CSAMiscProtocolHeaderInfoVA','CSANonImageHeaderInfoVB','CSAMiscProtocolHeaderInfoVB'};
% NB: might be necessary to add cases 'Private_0029_1110' and
% 'Private_0029_1210' for spectroscopic data (see spm_dicom_essentials.m)

for ccsa = 1:length(CSAlist)
    if isfield(hdr,CSAlist{ccsa})
        tdyhdr = tidy_CSA(hdr.(CSAlist{ccsa}));
        if isfield(tdyhdr,'MrPhoenixProtocol')
            % works for Siemens VB, VD & VE DICOM format:
            tdyhdr.MrPhoenixProtocol = read_ASCII(tdyhdr.MrPhoenixProtocol);
        elseif isfield(tdyhdr,'MrProtocol')
            % works for Siemens VA DICOM format:
            tdyhdr.MrProtocol = read_ASCII(tdyhdr.MrProtocol);
        end
        hdr.(CSAlist{ccsa}) = tdyhdr;
    end
end