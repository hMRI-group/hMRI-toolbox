%% Unittest for hMRI-toolbox versus VBQ toolbox
% =================
% How does it work?
% =================
% You will be asked for an output folder and your input data. 
% Afterwards the maps will be written into the given folder, calculated
% first with the old VBQ- and than with the hMRI-toolbox.
% Make sure that the toolbox settings are chosen similar.
% After the maps are calculated, brain masked differences will be
% calculated and deviations from 0 will be reported.
% Written by Tobias Leutritz

nrun = 1; % enter the number of runs here
[scriptpath,~,~] = fileparts(mfilename('fullpath'));
jobfile_VBQ = {strcat(scriptpath,filesep,'UNICORT_VBQ_job.m')};
jobfile_hMRI = {strcat(scriptpath,filesep,'UNICORT_job.m')};
jobs_VBQ = repmat(jobfile_VBQ, 1, nrun);
jobs_hMRI = repmat(jobfile_hMRI, 1, nrun);
inputs = cell(4, nrun);
for crun = 1:nrun
    inputs{1, crun} = cfg_getfile(Inf,'dir','Select output folder');
    inputs{2, crun} = cellstr(spm_select(Inf,'image','Select raw multiparameter data: MT images'));
    inputs{3, crun} = cellstr(spm_select(Inf,'image','Select raw multiparameter data: PD images'));
    inputs{4, crun} = cellstr(spm_select(Inf,'image','Select raw multiparameter data: T1 images'));
end
spm('defaults', 'FMRI');

spm_jobman('run', jobs_VBQ, inputs{:});

% override defaults to prevent MPMcalc to be deleted and ACPC realignment
% since the VBQ-toolbox doesn't realign to MNI-space
global hmri_def
hmri_def.cleanup = false;
hmri_def.qMRI_maps.ACPCrealign = 0;

spm_jobman('run', jobs_hMRI, inputs{:});

%% doing the unit-test afterwards
% output versions of the toolbox, SPM, Matlab, OS, Hardware (processor,
% RAM)

% define maps to be compared and read in for brain mask calculation
% #              1    2            3     4        5     6 
hMRI_map.name = {'R1' 'R1_masked' 'R2s' 'R2s_OLS' 'MT' 'A' ...
    'MT_outer_suppressed' 'MT_outer_suppressed'};
%     7                      8     

hMRI_map.name_prefix = {'' 'm' '' '' '' '' 'c1' 'c2'}; 
% #                     1   2  3  4  5   6   7   8  

VBQ_map.name = {'R1' 'R1' 'R2s' 'R2s_OLS' 'MT' 'A'};
% #              1    2    3     4        5     6 

VBQ_map.name_prefix = {'' 'mh' '' '' '' ''}; 
% #                     1   2  3  4  5   6 

map.name = {'R1' 'R1_UNICORT' 'R2s' 'R2s_OLS' 'MT' 'A'};

nmc = numel(hMRI_map.name);
omc = numel(VBQ_map.name);

for crun = 1:nrun

    % initialise cell variables
    vol_info_hMRI = cell(nmc,1);
    vol_info_VBQ = cell(omc,1);
    mnames_hMRI = cell(omc,1);
    mnames_VBQ = cell(omc,1);
    
    %% read in VBQ-toolbox results
    cd(char(inputs{1, crun}));
    % extract basename of files
    bn_key = '_A.nii';
    bnf = dir(strcat('*',bn_key));
    bn = strrep(bnf.name,bn_key,'');

    % loop over given number of maps
    for mn = 1:omc
        % create file names of maps
        mnames_VBQ{mn} = strcat(cd,filesep,VBQ_map.name_prefix{mn}, ...
            bn,'_',VBQ_map.name{mn},'.nii');
        % create SPM handles
        vol_info_VBQ{mn,1} = spm_vol(mnames_VBQ{mn});
    end
    
    %% read in hMRI-toolbox results
    cd(char(strcat(inputs{1, crun},filesep,'MPMCalc')));
    % extract basename of files
    bn_key = '_A.nii';
    bnf = dir(strcat('*',bn_key));
    bn = strrep(bnf.name,bn_key,'');

    % loop over given number of maps
    for mn = 1:nmc
        % create file names of maps
        mnames_hMRI{mn} = strcat(cd,filesep,hMRI_map.name_prefix{mn}, ...
            bn,'_',hMRI_map.name{mn},'.nii');
        % create SPM handles
        vol_info_hMRI{mn,1} = spm_vol(mnames_hMRI{mn});
    end
    
    % create brain mask from GM / WM segmentations
    brain_mask = vol_info_hMRI{1,1};
    brain_mask.fname = 'brain_mask.nii';
    brain_mask = spm_imcalc([vol_info_hMRI{7,1} vol_info_hMRI{8,1}], ...
        brain_mask,'(i1>0.9)+(i2>0.9)');
    
    % compare masked maps
    for mn = 1:omc
        map_name = map.name{mn};
        diff.(map_name) = vol_info_hMRI{1,1};
        diff.(map_name).fname = strcat(cd,filesep,'diff_',map_name,'.nii');
        diff.(map_name) = spm_imcalc([brain_mask vol_info_VBQ{mn,1} ...
            vol_info_hMRI{mn,1}],diff.(map_name),'(i1.*i2)-(i1.*i3)');
        diff_data.(map_name) = spm_read_vols(diff.(map_name));
        if isempty(nonzeros(diff_data.(map_name))) 
            delete(diff.(map_name).fname);
        else
            fprintf(1,'Differences in %s map found. Please check file %s!', ...
                map_name, diff.(map_name).fname)
        end
    end

end

