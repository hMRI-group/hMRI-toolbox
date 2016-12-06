function version = hMRI_get_version
% To retrieve the SHA1, author, date and message of the last commit of the
% current branch of the repository and return the information as a string.
% This script MUST be located in the root directory of the repository.
% If the Toolbox has been copied whitout version tracking, the version can
% only be retrieved if a version.txt file is already present in the root
% directory of the Toolbox.
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

% initialise the output variable:
version = [];

try
    % execute the git command to output the current version into text file:
    git log -1 > version.txt
catch MExc %#ok<*NASGU>
    % warning('Either or both MATLAB-git and git are not available on this machine.');
    fprintf(1,['\nWARNING:\nEither or both MATLAB-git and git are not available on this machine.\n'...
        'The current version of the hMRI-Toolbox cannot be retrieved.\n'...
        'See ''help hmri_get_version'' for details about dependencies.\n'...     
        'Searching for an existing version.txt file...\n']);
    fprintf(1,'\n%s\n', MExc.getReport);
end

% version.txt file should be the following:
version_fname = fullfile(repos_dir,'version.txt');

if ~exist(version_fname,'file')
    fprintf(1,'File %s does not exist.\nhMRI-Toolbox version unknown.\n\n', version_fname);
    version = 'Unknown version';
else
    % try to open the version.txt file:
    fid = fopen(version_fname,'r');
    if (fid~=-1)
        % read file content
        clin = fgets(fid);
        while (clin~=-1)
            version = [version clin];
            clin = fgets(fid);
        end
        fclose(fid);
        % fprintf(1,'Toolbox version:\n%s\n', version);
    else
        fprintf(1,'Cannot open file %s.\nhMRI-Toolbox version unknown.\n\n', version_fname);
        version = 'Unknown version';
    end
end

% back to the current working directory
cd(current_dir);

end
