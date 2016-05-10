function eb_json_metadata_example_01
% NOTES
% - make sure the file
%   'json_metadata\examples\f720-0008-00094-000094-01.nii' is present
%   (might need to make a copy of the original file).
% - change the paths below are correctly set...

% toolbox directory
cdir = pwd;
% addpath to the metadata set/get scripts
addpath(fullfile(cdir, '..'));

% initial file has no extended header
filename = fullfile(cdir,'f720-0008-00094-000094-01.nii');

% read unexisiting extended header:
hdr = get_extended_hdr(filename);
disp(hdr{1});

% define a structure with metadata to be included in the nifti extended
% header:
newhdr = struct('TR',8300,'TE',72,'FA',90,'an_array',[43.9 0.32 89;1.2 47 78],'AcquisitionDate','16-Mar-2016');
disp(newhdr);

% insert this structure as metadata extended header into the nifti file
set_extended_hdr(filename,newhdr);

% At this stage, the file has grown a tiny bit bigger and the
% JSON-formatted header can be read when opening the nifti file with a text
% editor.

% NOTE: in practice, the metadata should be first read
% (get_extended_hdr), then updated and reinjected into the nifti file,
% in order to avoid overwriting previous metadata of interest!

% Check this has definitely been done properly and we can retrieve the structure:
newhdr = get_extended_hdr(filename);
disp(newhdr{1});

% Rewind (partly) what has been done:
set_extended_hdr(filename,[]);
