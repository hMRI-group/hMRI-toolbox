function version = hMRI_get_version
% To retrieve the SHA1, author, date and message of the last commit of the
% current branch of the repository and return the information as a string.
% This script MUST be located in the root directory of the repository.
%
% DEPENDENCIES
% This script calls the git command using the MATLAB-git wrapper from
% https://github.com/manur/MATLAB-git.git. The latter allows you to use
% command line git instructions in Matlab (as long as Git is installed on
% your computer!). Make sure that the git.m script is in the Matlab path to
% execute this script.
%
% COMMAND LINE EQUIVALENT IN GIT BASH 
% (to output the information into version.txt)
% § git log -1 > version.txt
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - May 2016
% Cyclotron Research Centre, University of Liege

% retrieve the directory containing the local repository (i.e. the
% directory containing the current script)
repos_dir = fileparts(mfilename('fullpath'));
% the git command must be run from the repository directory so...
% keep track of where we are before running the git command
current_dir = pwd;
% move to the repository
cd(repos_dir);
% execute the git commant
git log -1 > version.txt
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

% back to the current working directory
cd(current_dir);

end
