function out = spm_run_dicom(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_dicom.m 6376 2015-03-12 15:15:57Z john $


if ~isempty(job.outdir{1})
    out_dir = job.outdir{1};
else
    out_dir = pwd;
end

if job.convopts.icedims
    root_dir = ['ice' job.root];
else
    root_dir = job.root;
end

% determine whether JSON metadata are stored and whether it is as extended
% header and/or separate JSON file.
json.extended = false;
json.separate = false;
if strfind(job.convopts.metaopts.mformat,'ext')
    if strcmp(job.convopts.format,'nii')
        json.extended = true; 
    end
end
if strfind(job.convopts.metaopts.mformat,'sep')
    json.separate = true; 
end

% determine whether full content of header is dumped into the JSON metadata
% structure or only the "essentials":
essentials = true;
if job.convopts.metaopts.mcontent
    essentials = false;
end

% determine the degree of confidentiality / anonymisation
json.anonym = job.convopts.metaopts.manonym;

hdr = spm_dicom_headers(char(job.data), essentials);

sel = true(size(hdr));
if ~isempty(job.protfilter) && ~strcmp(job.protfilter, '.*')
    psel   = cellfun(@(h)isfield(h, 'ProtocolName'), hdr);
    ssel   = ~psel & cellfun(@(h)isfield(h, 'SequenceName'), hdr);
    pnames = cell(size(hdr));
    pnames(psel) = cellfun(@(h)subsref(h, substruct('.','ProtocolName')), hdr(psel), 'UniformOutput', false);
    pnames(ssel) = cellfun(@(h)subsref(h, substruct('.','SequenceName')), hdr(ssel), 'UniformOutput', false);
    sel(psel|ssel) = ~cellfun(@isempty,regexp(pnames(psel|ssel), job.protfilter));
end
out = spm_dicom_convert(hdr(sel),'all',root_dir,job.convopts.format,out_dir, json);

