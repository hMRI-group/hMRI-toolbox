function hmri_unittest(varargin)
%% Unittest for hMRI-toolbox
% =================
% How does it work?
% =================
% First of all this script requires SPM to be up and runnning.
% You will be asked for input folders containing output data from the old 
% toolbox version (first folder) and the new version (second folder), 
% whereas you point to the folder containing the folders B1mapCalc,
% MPMCalc, Results etc. or the plain maps in case of VBQ-toolbox versions.
% Brain masked differences will be calculated and deviations from 0 will 
% be reported in the file comp_results.txt in the folder of the new version.
% Some thresholds can be adapted (lines 17-19).
%
% Written by Tobias Leutritz (30.06.2017)

bm_thr = 0.99; % threshold for brain masking
pu_thr = 98; % upper threshold for percentiles
pl_thr = 2; % lower threshold for percentiles

if nargin == 0
    % reading data directories, containing the map data
    ddirs = cfg_getfile(inf,'dir','Select pairwise output folders to compare results');
elseif mod(nargin,2) ~= 0
    warning('There have to be 2n input directories to compare!')
    ddirs = cfg_getfile(inf,'dir','Select pairwise output folders to compare results');
else
    np = numel(varargin);
    for n = 1:np
        ddirs{n,1} = varargin{1,n}; %#ok<AGROW>
    end
end

for dc = 1:numel(ddirs)/2
%% read in old toolbox results
cd(char(ddirs(dc*2-1,:)));
new_ver = false; % check variable to determine if hMRI-Toolbox was used
if isdir('Results')
    try
        cd('MPMCalc');
        new_ver = true;
    catch
        fprintf(1,['No MPMCalc folder existing - please rerun the ' ...
            'calculation with hmri_defauls.cleanup = 0!\n']);
        return;
    end
end
% extract basename of files
bn_key = '_MT.nii';
bnf = dir(strcat('*',bn_key));
bn = strrep(bnf.name,bn_key,'');

% all NIFTI files
all_files = dir('*.nii');
all_files = {all_files.name};
nmpmc = numel(all_files);
if new_ver
    cd('..'); cd('B1mapCalc');
    all_files2 = dir('*.nii');
    all_files2 = {all_files2.name};
    all_files = [all_files, all_files2]; %#ok<AGROW>
    cd('..');
end
    
nm = numel(all_files);

% initialise cell variables
vol_info_new = cell(nm + 3,1); % +3 for c1, c2, c3 tissue segments
vol_info_old = cell(nm,1);
mnames_new = cell(nm + 3,1);
mnames_old = cell(nm,1);
mc = 0; % counter for missing maps

for n = 1:nm
    % create filename with full path
    if new_ver
        if n <= nmpmc
            mnames_old{n,1} = fullfile(cd,strcat('MPMCalc',filesep,all_files(n)));
        else
            mnames_old{n,1} = fullfile(cd,strcat('B1mapCalc',filesep,all_files(n)));
        end
    else
        mnames_old{n,1} = fullfile(cd,all_files(n));
    end
    % create SPM handles
    vol_info_old{n,1} = spm_vol(char(mnames_old{n,1}));
end

%% read in new toolbox results
cd(char(ddirs(2*dc,:)));
output_file = fopen('comp_results.txt','w');
spm_version = spm('Version');
mtl_version = version;
sys_version = computer;
fprintf(output_file,['Unittest for hMRI-Toolbox version %s computed with %s' ... 
    ' under MATLAB %s on a %s computer system.\n\n Following directories contain'...
    ' the data being compared:\n %s \n %s\n\n'], ...
     hmri_get_version, spm_version, mtl_version, sys_version, ...
     char(ddirs(dc*2-1,:)), char(ddirs(2*dc,:)));

cd(char(strcat(ddirs(2*dc,:),filesep,'MPMCalc')));

% loop over given number of maps
for mn = 1:nm
    % create file names of maps
    mnames_new{mn,1} = fullfile(cd,all_files(mn));
    % create SPM handles & check for VBQ-Toolbox file naming & correct that
    if cell2mat(strfind(mnames_new{mn,1},'_A.nii'))
        try
           vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        catch
           mnames_new{mn,1} = fullfile(cd,strrep(all_files(mn),'_A.nii','_PD.nii'));
           vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        end
    elseif cell2mat(strfind(mnames_new{mn,1},'_MTforA.nii'))
        try
           vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        catch
           mnames_new{mn,1} = fullfile(cd,strrep(all_files(mn),'_MTforA.nii','_MT_outer_suppressed.nii'));
           vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        end
    elseif (~isempty(cell2mat(strfind(mnames_new{mn,1},'_R1.nii'))) && ...
            ~isempty(cell2mat(strfind(mnames_new{mn,1},'B1_s'))))
        try
           vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        catch
            oldd = cd;
            cd(char(strcat(ddirs(2*dc,:),filesep,'B1mapCalc')));
            mnames_new{mn,1} = fullfile(cd,strrep(all_files(mn),'B1_s','B1u_s'));
            vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
            cd(oldd);
        end
    elseif (~isempty(cell2mat(strfind(mnames_new{mn,1},'_R1.nii'))) && ...
            ~isempty(cell2mat(strfind(mnames_new{mn,1},'hs'))))
        try
           vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        catch
            temp_name = strrep(all_files(mn),'h','');
            temp_name = strrep(temp_name,'R1','R1_masked');
            mnames_new{mn,1} = fullfile(cd,temp_name);
            vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        end
    else
        try
            vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
        catch 
           try % go to B1mapCalc folder and look there...
               oldd = cd;
               cd(char(strcat(ddirs(2*dc,:),filesep,'B1mapCalc')));
               mnames_new{mn,1} = fullfile(cd,all_files(mn));
               vol_info_new{mn,1} = spm_vol(char(mnames_new{mn,1}));
               cd(oldd);
           catch
               [~, map_name, ~] = fileparts(char(mnames_new{mn,1}));
                fprintf(output_file,'Not found among new hMRI-Toolbox results: %s\n', ...
                    strcat(map_name,'.nii'));
                mc = mc + 1;
                missing(mc) = mn; %#ok<SAGROW>
                cd(oldd);
           end
        end  
    end
end

%% create brain mask from GM, WM and CSF segmentations
% read in tissue segmentations
cm1 = nm + 1;
mnames_new{cm1,1} = fullfile(cd,strcat('c1',bn,'_MT_outer_suppressed.nii'));
vol_info_new{cm1,1} = spm_vol(mnames_new{cm1,1});
cm2 = nm + 2;
mnames_new{cm2,1} = fullfile(cd,strcat('c2',bn,'_MT_outer_suppressed.nii'));
vol_info_new{cm2,1} = spm_vol(mnames_new{cm2,1});
cm3 = nm + 3;
mnames_new{cm3,1} = fullfile(cd,strcat('c3',bn,'_MT_outer_suppressed.nii'));
vol_info_new{cm3,1} = spm_vol(mnames_new{cm3,1});

% prepare and write brain mask
cd(char(ddirs(2*dc,:)));
brain_mask = vol_info_new{nm + 1,1};
brain_mask.descrip = 'Brain Mask';
brain_mask.fname = strcat(cd,filesep,'brain_mask.nii');
brain_mask = spm_imcalc([vol_info_new{cm1,1} vol_info_new{cm2,1} ...
    vol_info_new{cm3,1}], brain_mask, ...
    sprintf('((i1>%f)+(i2>%f)+(i3>%f))>0',bm_thr, bm_thr, bm_thr));
masks = {'GM' 'WM' 'CSF'};
for n = 1:3
    mask_name = char(masks(n));
    cm  = eval(sprintf('cm%i',n));
    mask.(mask_name) = vol_info_new{cm,1};
    mask.(mask_name).descrip = strcat(mask_name,' Mask');
    mask.(mask_name).fname = strcat(cd,filesep,mask_name,'.nii');
    mask.(mask_name) = spm_imcalc(vol_info_new{cm,1},  ...
        mask.(mask_name), sprintf('i1>%f',bm_thr));
end

%% compare masked maps
for mn = 1:nm
    % skip missing files
    if isempty(vol_info_new{mn,1})
        continue;
    end
    [pth, map_name, ~] = fileparts(char(mnames_new{mn,1}));
    short_name = strcat(strrep(map_name,bn,'[...]'),'.nii');
    diff = vol_info_new{mn,1};    
    diff.descrip = 'Absolute Difference';    
	diff.fname = strcat(cd,filesep,'diff_',short_name);
    % skip brain masking for B1map related files
    if strfind(pth,'B1mapCalc')
        diff = spm_imcalc([vol_info_old{mn,1} vol_info_new{mn,1}], ...
            diff,'i2 - i1');
        bm = false;
        else
        diff = spm_imcalc([brain_mask vol_info_old{mn,1} ...
            vol_info_new{mn,1}],diff,'(i1.*i3)-(i1.*i2)');
        bm = true;
    end
    % read data and look for nonzeros
    diff_data = spm_read_vols(diff);
    temp_data = nonzeros(diff_data);
    temp_data(isnan(temp_data)) = [];
    % output results
    if isempty(nonzeros(temp_data))
        fprintf(output_file,'No differences found in "%s" (%s).\n', ...
            vol_info_new{mn,1}.descrip, short_name);
        delete(diff.fname);
    else
        % calculate relative difference
        div = vol_info_new{mn,1};    
        div.descrip = 'Relative Difference [p.u.]';    
        div.fname = strcat(cd,filesep,'div_',short_name);
        div = spm_imcalc([diff vol_info_old{mn,1}],div,'i1./i2*100');
        div_data = spm_read_vols(div);
        temp_data = nonzeros(div_data);
        temp_data(isnan(temp_data)) = [];
        minp.bm = prctile(temp_data(:),pl_thr); % min(temp_data(:));
        maxp.bm = prctile(temp_data(:),pu_thr); % max(temp_data(:));
        if ~bm
            fprintf(output_file,['Differences found in "%s" (%s) in a range of'...
                ' [%.2f%%; %.2f%%]. Please check files "%s" and "%s"!\n'], ...
                vol_info_new{mn,1}.descrip, short_name, minp.bm, maxp.bm, ...
                strcat('diff_',map_name,'.nii'), strcat('div_',map_name,'.nii'));
        else
            for n = 1:3
                mask_name = char(masks(n));
                temp_vol = div;
                temp_vol.fname = strcat(cd,filesep,'temp.nii');
                temp_vol = spm_imcalc([div mask.(mask_name)],temp_vol,'i1.*i2');
                temp_data = spm_read_vols(temp_vol);
                temp_data = nonzeros(temp_data);
                temp_data(isnan(temp_data)) = [];
                minp.(mask_name) = prctile(temp_data(:),pl_thr); % min(temp_data(:));
                maxp.(mask_name) = prctile(temp_data(:),pu_thr); % max(temp_data(:));
            end
            fprintf(output_file,['Differences found in "%s" (%s) in a total range of'...
                ' [%.2f%%; %.2f%%] (GM: [%.2f%%; %.2f%%], WM: [%.2f%%; %.2f%%]' ...
                ', CSF: [%.2f%%; %.2f%%]). Please check files "%s" and "%s"!\n'], ...
                vol_info_new{mn,1}.descrip, short_name, minp.bm, maxp.bm, ...
                minp.GM, maxp.GM, minp.WM, maxp.WM, minp.CSF, maxp.CSF, ...
                strcat('diff_',map_name,'.nii'), strcat('div_',map_name,'.nii'));
        end
    end
end

fclose(output_file);
delete('temp.nii');

end