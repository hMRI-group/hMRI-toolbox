function out = spm_run_dicom(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_dicom.m 7201 2017-11-08 11:13:25Z guillaume $
% modified by Evelyne Balteau for hMRI-toolbox compatibility

% to make sure only hMRI implementation is used
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

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
% structure or only the "essentials" (default is false, i.e. full content):
essentials = false;
if isfield(job.convopts.metaopts,'mcontent')
    if ~job.convopts.metaopts,mcontent
        essentials = true;
    end
end

% determine the degree of confidentiality / anonymisation (attempt!)
json.anonym = 'basic';
if isfield(job.convopts.metaopts,'manonym')
    json.anonym = job.convopts.metaopts.manonym;
end

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

% to resolve conflicts with SPM version, add path to hMRI so it goes back
% to top of the list if it was not the case...
addpath(fileparts(fileparts(mfilename('fullpath'))));
out = spm_dicom_convert(hdr(sel),'all',root_dir,job.convopts.format,out_dir, json);

