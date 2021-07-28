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


function tests = test_find_field(testCase)
  inStruct = testCase.TestData.testJSON;

  % Non-recursive search
  [res, status] = hmri_metadata.find_field(inStruct, 'topLevelScalar', true);
  assert(status);
  assert(res == 2);

  [res, status] = hmri_metadata.find_field(inStruct, 'toplevelscalar', false);
  assert(status);
  assert(res == 2);

  [res, status] = hmri_metadata.find_field(inStruct, 'topLevelString', true);
  assert(status);
  assert(strcmp(res, 'abc'));

  [res, status] = hmri_metadata.find_field(inStruct, 'toplevelstring', true);
  assert(~status);
  assert(isempty(res));

  [res, status] = hmri_metadata.find_field(inStruct, ...
                                           sprintf('  topLevelString\t'), ...
                                           true);
  assert(status);
  assert(strcmp(res, 'abc'));


  [res, status] = hmri_metadata.find_field(inStruct, 'topLevelList', true);
  assert(status);
  assert(all(res == [1; 2; 3]));

  [res, status] = hmri_metadata.find_field(inStruct, 'topLevelLongList', true);
  assert(status);
  assert(all(res == [11; 12; 13; 21; 22; 23]));

  [res, status] = hmri_metadata.find_field(inStruct, 'topLevelStruct', true);
  assert(status);
  assert(res.subLevelScalar == 3.14);
  assert(strcmp(res.subLevelString, 'def'));

  % Recursive search
  res = hmri_metadata.find_field(inStruct, 'topLevelNested', true);
  res
  [res, status] = hmri_metadata.find_field(inStruct, 'item1', true);
  assert(status);
  assert(strcmp(res.i1j1, 'val_i1j1'));


  [res, status] = hmri_metadata.find_field(inStruct, 'i1j3', true);
  assert(status);
  assert(strcmp(res, 'val_i1j3'));
end
