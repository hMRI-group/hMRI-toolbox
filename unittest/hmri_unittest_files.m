function hmri_unittest_files(varargin)
%% Unittest for hMRI-toolbox
% =================
% How does it work?
% =================
% First of all this script requires SPM to be up and runnning.
% You will be asked for input files containing output data from the old 
% toolbox version (first set) and the new version (second set).
% No differenciation in filenames will be done, i.e. the order of input
% files has to be the same for both sets of data.
% Differences will be calculated and deviations from 0 will 
% be reported in the file comp_results.txt in the folder of the new version.
% Some thresholds can be adapted (lines 15-16).
%
% Written by Tobias Leutritz (30.06.2017)

pu_thr = 98; % upper threshold for percentiles
pl_thr = 2; % lower threshold for percentiles

if nargin == 0
    % reading data directories, containing the map data
    mnames_old = cfg_getfile(inf,'image','Select output files of the old version');
    mnames_new = cfg_getfile(inf,'image','Select output files of the new version');
elseif mod(nargin,2) ~= 0
    warning('There have to be 2n input files to compare!');
    mnames_old = cfg_getfile(inf,'image','Select output files of the old version');
    mnames_new = cfg_getfile(inf,'image','Select output files of the new version');
end


while numel(mnames_old) ~= numel(mnames_new)
    warning(['Mismatch in number of given input. ' ...
        'Please input same amount of files in the following dialogs!']);
    mnames_old = cfg_getfile(inf,'image','Select output files of the old version');
    mnames_new = cfg_getfile(inf,'image','Select output files of the new version');
end

nm = numel(mnames_old);

%% read in results

% initialise cell variables
vol_info_new = cell(nm,1);
vol_info_old = cell(nm,1);

for n = 1:nm
    % create SPM handles
    vol_info_old{n,1} = spm_vol(char(mnames_old{n}));
    vol_info_new{n,1} = spm_vol(char(mnames_new{n}));
end

%% write info file
[pth1, ~, ~] = fileparts(mnames_old{1});
[pth2, ~, ~] = fileparts(mnames_new{1});
cd(pth2);
output_file = fopen('comp_results.txt','w');
spm_version = spm('Version');
mtl_version = version;
sys_version = computer;
fprintf(output_file,['Unittest for hMRI-Toolbox version %s computed with %s' ... 
    ' under MATLAB %s on a %s computer system.\n\n Following directories contain'...
    ' the data being compared:\n %s \n %s\n\n'], ...
     hmri_get_version, spm_version, mtl_version, sys_version, pth1, pth2);

%% compare maps
for mn = 1:nm
    [~, map_name, ~] = fileparts(char(mnames_new{mn}));
    diff = vol_info_new{mn,1};    
    diff.descrip = 'Absolute Difference';    
	diff.fname = strcat(cd,filesep,'diff_',map_name,'.nii');
    diff = spm_imcalc([vol_info_old{mn,1} vol_info_new{mn,1}], diff, 'i2-i1');
    % read data and look for nonzeros
    diff_data = spm_read_vols(diff);
    temp_data = nonzeros(diff_data);
    temp_data(isnan(temp_data)) = [];
    % output results
    if isempty(nonzeros(temp_data))
        fprintf(output_file,'No differences found in "%s" (%s).\n', ...
            vol_info_new{mn,1}.descrip, map_name);
        delete(diff.fname);
    else
        % calculate relative difference
        div = diff;
        div.descrip = 'Relative Difference [p.u.]';    
        div.fname = strcat(cd,filesep,'div_',map_name,'.nii');
        div = spm_imcalc([diff vol_info_old{mn,1}],div,'i1./i2*100');
        div_data = spm_read_vols(div);
        temp_data = nonzeros(div_data);
        temp_data(isnan(temp_data)) = [];
        minp = prctile(temp_data(:),pl_thr); % min(temp_data(:));
        maxp = prctile(temp_data(:),pu_thr); % max(temp_data(:));
        fprintf(output_file,['Differences found in "%s" (%s) in a range of'...
                ' [%.2f%%; %.2f%%]. Please check files "%s" and "%s"!\n'], ...
                vol_info_new{mn,1}.descrip, map_name, minp, maxp, ...
                strcat('diff_',map_name,'.nii'), strcat('div_',map_name,'.nii'));
    end
end

fclose(output_file);
delete('temp.nii');

end