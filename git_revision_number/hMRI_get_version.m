function version = hMRI_get_version
% Read the information stored in version.txt (located in the root directory
% of the repository) and return it as a string. Used when creating/updating
% the extended header with the version of the software at the time of the
% creation/update.
%
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - May 2016
% Cyclotron Research Centre, University of Liege

% retrieve the directory containing the local repository (i.e. the
% directory containing the current script)
repos_dir = fileparts(mfilename('fullpath'));
% open the version.txt file
version_fname = fullfile(repos_dir,'version.txt');
fid = fopen(version_fname,'r');

version = [];

if (fid~=-1)
    % read file content
    clin = fgets(fid);
    while (clin~=-1)
        version = [version clin];
        clin = fgets(fid);
    end
    fclose(fid);
else
    warning('Cannot open file %s.', version_fname);
end