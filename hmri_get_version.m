function hmri_version = hmri_get_version
% To retrieve the current hMRI-Toolbox version number from version.txt.
% Additional information (SHA1, author, date and message of the last commit
% of the current branch of the repository) can be retrieved using git from
% Matlab (optional, see DEPENDENCIES below).
% This script MUST be located in the root directory of the repository.
% If the Toolbox has been copied whitout version tracking, the version can
% only be retrieved if a version.txt file is already present in the root
% directory of the Toolbox. No additional commit number will be retrieved.
%
% DEPENDENCIES (NOT MANDATORY)
% This script calls the git command using the MATLAB-git wrapper from
% https://github.com/manur/MATLAB-git.git. The latter allows you to use
% command line git instructions in Matlab (as long as Git is installed on
% your computer! - https://git-scm.com/download). Make sure that the git.m
% script is in the Matlab path to execute this script. 
%
% COMMAND LINE EQUIVALENT IN GIT BASH 
% (to output the information into lastcommit.txt)
% § git log -1 > lastcommit.txt
%
% ALTERNATIVELY, see version/hmri_get_version_readme.pdf for alternative
% way to update the version.txt file using post-commit and post-merge
% scripts.
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

% initialise the output variable:
hmri_version = [];

try
    % execute the git command to output the current version into text file:
    git log -1 > lastcommit.txt
    lastcommit = true;
catch MExc %#ok<*NASGU>
    lastcommit = false;
    % fprintf(1,['\nWARNING:\nEither or both MATLAB-git and git are not available on this machine.\n'...
    %    'The current version of the hMRI-Toolbox cannot be retrieved.\n'...
    %    'See ''help hmri_get_version'' for details about dependencies.\n'...     
    %    'Searching for an existing version.txt file...\n']);
    % fprintf(1,'\n%s\n', MExc.getReport);
end

% version.txt file should be the following:
version_fname = fullfile(repos_dir,'version.txt');
% lastcommit.txt file should be the following:
lastcommit_fname = fullfile(repos_dir,'lastcommit.txt');

if ~exist(version_fname,'file')
    fprintf(1,'File %s does not exist.\nhMRI-Toolbox version unknown.\n\n', version_fname);
    hmri_version = 'Unknown hMRI version. File version.txt does not exist.';
else
    % try to open the version.txt file:
    fid = fopen(version_fname,'r');
    if (fid~=-1)
        % init version description
        hmri_version = 'hMRI ';
        % read file content
        clin = fgets(fid);
        while (clin~=-1)
            hmri_version = [hmri_version clin]; %#ok<AGROW>
            clin = fgets(fid);
        end
        fclose(fid);
        % fprintf(1,'Toolbox version:\n%s\n', version);
    else
        fprintf(1,'Cannot open file %s.\nhMRI-Toolbox version unknown.\n\n', version_fname);
        hmri_version = 'Unknown hMRI version. Cannot open file version.txt.';
    end
end

% add last commit if available 
if lastcommit
    % try to open the version.txt file:
    fid = fopen(lastcommit_fname,'r');
    if (fid~=-1)
        % read file content
        clin = fgets(fid);
        while (clin~=-1)
            hmri_version = [hmri_version clin]; %#ok<AGROW>
            clin = fgets(fid);
        end
        fclose(fid);
        % fprintf(1,'Toolbox version:\n%s\n', version);
    end
    delete(lastcommit_fname);
end

% add Matlab version for full version tracking
hmri_version = sprintf('Matlab %s \n%s', version, hmri_version);

% back to the current working directory
cd(current_dir);

end
