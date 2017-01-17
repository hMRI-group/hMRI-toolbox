function hdrout = spm_dicom_anonymise(hdr, opts)
% USAGE: hdrout = spm_dicom_anonymise(hdr, opts)
% hdr is a DICOM header read with spm_dicom_header
% opts are anonymisation options:
%       opts.anonym = 'full': no patient information is kept at all
%                     'basic': patient ID (presumably not containing his
%                     name), age (years at the time of the data
%                     acquisition), sex, size and weight are kept, patient
%                     name, date of birth and DICOM filename (often
%                     containing the patient name) are removed.
%=========================================================================%
% Evelyne Balteau - Cyclotron Research Centre - May 2016
%=========================================================================%

if nargin<2
    opts.anonym = 'basic';
end
if ~iscell(hdr)
    hdrout{1} = hdr;
else
    hdrout = hdr;
end

for i=1:length(hdrout)
    if isfield(hdrout{i},'PatientName')
        hdrout{i}.PatientName = 'Anonymous';
    end
    if isfield(hdrout{i},'PatientBirthDate')
        t1 = hdrout{i}.PatientBirthDate;
        t2 = hdrout{i}.StudyDate;
        hdrout{i}.PatientAge = round((t2-t1)*10/365.25)/10;
        hdrout{i} = rmfield(hdrout{i},'PatientBirthDate');
    end
    if isfield(hdrout{i},'Filename')
        hdrout{i} = rmfield(hdrout{i},'Filename');
    end
    
    if strcmp(opts.anonym, 'full')
        try hdrout{i} = rmfield(hdrout{i},'PatientID');end %#ok<*TRYNC>
        try hdrout{i} = rmfield(hdrout{i},'PatientSex');end
        try hdrout{i} = rmfield(hdrout{i},'PatientAge');end
        try hdrout{i} = rmfield(hdrout{i},'PatientSize');end
        try hdrout{i} = rmfield(hdrout{i},'PatientWeight');end
        try hdrout{i}.CSASeriesHeaderInfo = rmfield(hdrout{i}.CSASeriesHeaderInfo,'UsedPatientWeight');end
    end
end