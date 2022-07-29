% Retrieves data from the cloudshare hosted at Leipzig
% Currently downloads and overwrites local data if it exists.

function hmri_get_ut_data

% Location of test data for hMRI toolbox (Callaghan et al. DIB 2019):
url_loc = 'https://datashare.mpcdf.mpg.de/s/74KYwDihKWazqbg/download?path=%2Fhmri_sample_dataset_with_maps%2Fgre_field_mapping_1acq_rl_0005';
field_map_1 = '&files=anon_s2018-02-28_18-26-185100-00001-00001-1.nii';

% Location to write data to locally:
ut_data_dir = [fileparts(which('hmri_test_utils')) filesep 'example_data'];
[~,~] = mkdir(ut_data_dir);

% Download data for unit testing (overwriting any previous file)
%urlwrite([url_loc field_map_1], [ut_data_dir filesep 'field_map_1.nii']);
x=webread([url_loc field_map_1]);
fid=fopen([ut_data_dir filesep 'field_map_1.nii'],'w');
fwrite(fid,x);
fclose(fid);

end