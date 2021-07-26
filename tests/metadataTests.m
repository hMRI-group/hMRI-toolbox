function tests = metadataTests
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  basePath = fileparts(mfilename('fullpath'));
  bidsDataset = fullfile(basePath, 'data', 'bidsified');
  bidsData = cellstr(spm_select('FPListRec', bidsDataset, '.*\.nii'));

  fid = fopen(fullfile(basePath, 'metafields_bids.txt'));
  meta = textscan(fid, '%s');
  MetaFields = meta{1};
  fclose(fid);

  sourceDataset = fullfile(basePath, 'data', ...
                           'hmri_sample_dataset_with_maps');
  sourceDataset = spm_select('FPListRec', sourceDataset, '.*\.nii');
  sourceDataset = cellstr(sourceDataset);

  testCase.TestData.MetaFields = MetaFields;
  testCase.TestData.bidsData = bidsData;
  testCase.TestData.sourceDataset = sourceDataset;
end

function tests = testMetafields(testCase)

  for iFile = 1:size(testCase.TestData.bidsData, 1)
    fname = testCase.TestData.bidsData{iFile};
    [path, basename, ~] = fileparts(fname);
    fprintf('%s:\n', basename);
    bids_header = get_metadata(fname);
    bids_header = bids_header{1};

    source_name = get_metadata_val(bids_header, 'OriginalFile');
    assert(~isempty(source_name), 'Unable to get original file name');
    fprintf('\tOriginal: %s\n', source_name);
    source_index = find(endsWith(testCase.TestData.sourceDataset, ...
                                 source_name));

    assert(size(source_index, 1) > 0, ...
           'Can''t find original file');
    source_file = testCase.TestData.sourceDataset{source_index(1)};
    source_header = get_metadata(source_file);
    source_header = source_header{1};

    for i = 1:size(testCase.TestData.MetaFields, 1)
      meta = testCase.TestData.MetaFields{i, 1};
      source_val = get_metadata_val(source_header, meta);
      bids_val = get_metadata_val(bids_header, meta);

      % bypass of field strenght calculation
      if strcmp(meta, 'FieldStrength')
        assertEqual(bids_val, 3);
        continue;
      end

      % bypass not filled FA in TB1EPI
      if strcmp(meta, 'FlipAngle') && source_val == 0
        continue;
      end

      if ischar(source_val)
        assert(strcmp(deblank(source_val), deblank(bids_val)));
      else
        assert(all(abs(source_val - bids_val) < 1e-8));
      end
    end
  end
end
