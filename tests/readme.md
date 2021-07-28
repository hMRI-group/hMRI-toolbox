### Tests

Tests requires test data downloaded into `tests/data` folder:

- `hmri_sample_dataset_with_maps` containing full qMRI dataset with rwa and pre-processed data
- `bidsified` a bidsified version of qMRI dataset


#### Metadata retrieval tests

Tests in `metadataTests.m` runs over all data files in bidsified dataset, 
retrives the metadata fields from `metafields_bids.txt`, and compares them
with values retrieved from raw dataset.

Only values that are not empty in raw dataset are concidered.


```matlab
res = runtests('metadataTests.m')
```
