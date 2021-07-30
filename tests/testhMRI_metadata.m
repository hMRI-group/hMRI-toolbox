function tests = testhMRI_metadata
  tests = functiontests(localfunctions);
end


function setupOnce(testCase)

  basePath = fileparts(mfilename('fullpath'));
  dataPath = fullfile(basePath, 'data');
  bidsDataset = fullfile(dataPath, 'hmri_sample_dataset_bids');
  bidsData = cellstr(spm_select('FPListRec', bidsDataset, '.*\.nii'));

  % Adding generic paths
  %
  % matlabbatch for running spm_select
  addpath(fullfile(spm('dir'), 'matlabbatch'));
  % get_metadata series of functions
  addpath(fullfile(basePath, '..', 'spm12', 'metadata'));

  sourceDataset = fullfile(dataPath, 'hmri_sample_dataset_with_maps');
  sourceDataset = spm_select('FPListRec', sourceDataset, '.*\.nii');
  sourceDataset = cellstr(sourceDataset);

  testCase.TestData.bidsData = bidsData;
  testCase.TestData.sourceDataset = sourceDataset;

  testCase.TestData.testJSON = spm_jsonread(...
                               fullfile(dataPath, 'testJSON.json'));
end


function tests = test_get_val_raw(testCase)
  inStruct = testCase.TestData.testJSON;

  % Non-recursive search
  [res, status] = hmri_metadata.get_val_raw(inStruct, 'topLevelScalar', true);
  assert(status);
  assert(res == 2);

  [res, status] = hmri_metadata.get_val_raw(inStruct, 'toplevelscalar', false);
  assert(status);
  assert(res == 2);

  [res, status] = hmri_metadata.get_val_raw(inStruct, 'topLevelString', true);
  assert(status);
  assert(strcmp(res, 'abc'));

  [res, status] = hmri_metadata.get_val_raw(inStruct, 'toplevelstring', true);
  assert(~status);
  assert(isempty(res));

  [res, status] = hmri_metadata.get_val_raw(inStruct, ...
                                           sprintf('  topLevelString\t'), ...
                                           true);
  assert(status);
  assert(strcmp(res, 'abc'));


  [res, status] = hmri_metadata.get_val_raw(inStruct, 'topLevelList', true);
  assert(status);
  assert(all(res == [1; 2; 3]));

  [res, status] = hmri_metadata.get_val_raw(inStruct, 'topLevelLongList', true);
  assert(status);
  assert(all(res == [11; 12; 13; 21; 22; 23]));

  [res, status] = hmri_metadata.get_val_raw(inStruct, 'topLevelStruct', true);
  assert(status);
  assert(res.subLevelScalar == 3.14);
  assert(strcmp(res.subLevelString, 'def'));

  % Recursive search
  res = hmri_metadata.get_val_raw(inStruct, 'topLevelNested', true);
  [res, status] = hmri_metadata.get_val_raw(inStruct, 'item1', true);
  assert(status);
  assert(strcmp(res.i1j1, 'val_i1j1'));


  [res, status] = hmri_metadata.get_val_raw(inStruct, 'i1j3', true);
  assert(status);
  assert(strcmp(res, 'val_i1j3'));
end


function tests = test_find_bids_entity(testCase)
  fname = 'sub-aaa_ses-bbb_acq-ddd_sufXYZ.nii.gz';

  res = hmri_metadata.find_bids_entity(fname, 'sub-');
  assertEqual(testCase, res, 'aaa');

  res = hmri_metadata.find_bids_entity(fname, 'ses-');
  assertEqual(testCase, res, 'bbb');

  res = hmri_metadata.find_bids_entity(fname, 'suffix');
  assertEqual(testCase, res, 'sufXYZ');

  res = hmri_metadata.find_bids_entity(fname, 'extension');
  assertEqual(testCase, res, '.nii.gz');

  res = hmri_metadata.find_bids_entity(fname, 'xyz-');
  assertEqual(testCase, res, []);

end
