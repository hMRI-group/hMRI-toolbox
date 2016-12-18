function hdr = spm_dicom_anonymise(hdr, opts)
% USAGE: hdr = spm_dicom_anonymise(hdr, opts)
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

if isfield(hdr,'PatientName')
    hdr.PatientName = 'anonymous';
end
if isfield(hdr,'PatientBirthDate')
    t1 = datenum(hdr.PatientBirthDate,'yyyymmdd');
    t2 = datenum(hdr.StudyDate,'yyyymmdd');
    hdr.PatientAge = round((t2-t1)*10/365.25)/10;
    hdr = rmfield(hdr,'PatientBirthDate');
end
if isfield(hdr,'Filename')
    hdr.Filename = 'AnonymousFileName';
end

% if strcmp(opts.anonym, 'full')
%     if isfield(hdr,'PatientID');hdr = rmfield(hdr,'PatientID');end
%     if isfield(hdr,'PatientID');hdr = rmfield(hdr,'PatientID');end
%     if isfield(hdr,'PatientID');hdr = rmfield(hdr,'PatientID');end
%     if isfield(hdr,'PatientID');hdr = rmfield(hdr,'PatientID');end
%     if isfield(hdr,'PatientID');hdr = rmfield(hdr,'PatientID');end
% end
