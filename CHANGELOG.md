# hMRI-toolbox changelog

This changelog documents all notable changes to the hMRI-toolbox.

Most recent version numbers *should* follow the [Semantic Versioning](https://semver.org/spec/v2.0.0.html) principles (e.g. bug fixes: x.x.1 > x.x.2, new feature with backward compatibility: x.2.x > x.3.0, major release affecting the way data are handled and processed: 1.x.x > 2.0.0).

## [v0.4.0]

### Added
- weighted least squares R2* fitting
- cleaner input/output functions for map creation
- several unit tests
- nulling of bright voxels during unified segmentation to improve segmentation of anonymised data

## [v0.3.0]

### Added
- QUIQI
- quality control tool

## [v0.2.6]

### Fixed
- Make cell processing code run for all paths of `get_metadata_val` so that TR, TE and FA can be read from the nifti description field if needed

## [v0.2.5]

- Provides option to remove the small angle approximation when calculating R1 and PD
- Bugfix for reverse phase encoding direction 3D-EPI on Siemens scanners
- Experimental BIDS support

## [v0.2.4]

- added back accidentally removed B1 mapping options

## [v0.2.3]

- allow calculation of B1 maps independently of the rest of the quantitative maps

## [v0.2.2]

- fix problem with broken JSON serialisation
- open files read only where possible

## [v0.2.1]

fix for issue #6

## [v0.2.0]

### Fixed

Many small bugs - refer to commits for details...

### Added

- **DICOM import for Philips data**:           
    Properly accounts for the rescaling factor for quantitative analysis
    (see Chenevert *et al.* 2014).

- **Tracking Matlab version number** in the JSON metadata.

### Changed

- **MT saturation map file name** from `*_MT` to `*_MTsat`.

## [v0.1.2-beta]

### Fixed

Many small bugs - refer to commits for details...

### Added

- **Single echo VFA approach**:
Implemented with a default minimum number of echoes for R2* calculation of 4.
The number of echoes used is logged and a general warning informs the user about data interpretation when R2* map derived from only a few echoes.

- **Metadata & processing parameters**:
Bunch of modifications to improve the readability and completeness of the metadata, for each output image as well as for the processing parameters.

- **UNICORT-generated R1 and B1 to calculate PD and MT maps**:
This option has not been thoroughly tested.
Made available under "ADVANCED USERS ONLY": the option to use the R1-UNICORT-derived B1 map for B1 transmit bias correction in PD and MT maps can only be enabled by an advanced user implementing customized defaults.

- **Example files**: For defaults customization and toolbox configuration examples.

- **New option to disable the coregistration steps**:
Coregistration shouldn't be disabled but it can be convenient in specific cases (simulated data, phantom data).
Option made available under "ADVANCED USERS ONLY", i.e. can only be modified by an advanced user implementing customized defaults.
When enabled (defaults), all input images to hmri_create_MTProt.m (transmit and receive fields, T1w and MTw images) are coregistered to the PDw average (or TE=0 fit) image (see hmri_create_MTProt.m).

### Changed

- **Options for RF sensitivity bias field correction**:
The available options are now `None`, `Unified Segmentation`(default), `Single` and `Per contrast`.
Beware that the previously implemented `None`option corresponds to the current `Unified Segmentation`!!!
Batch files saved with the previous version are not fully compatible with the current version and must be adapted.

## [v0.1.1-beta]

First public beta-version of the hMRI-toolbox.

Modular structure including the following modules:
- `Configure toolbox`
- `DICOM Import`
- `Auto-Reorient`
- `Create hMRI maps`
- `Process hMRI maps`

## [v0.2.0] (released 2018-12-20)

### Fixed

Many small bugs - refer to commits for details...

### Added

- **DICOM import for Philips data**:           
    Properly accounts for the rescaling factor for quantitative analysis
    (see Chenevert *et al.* 2014).

- **Tracking Matlab version number** in the JSON metadata.

### Changed

- **MT saturation map file name** from `*_MT` to `*_MTsat`.


## [v0.1.2-beta2] (released 2018-07-30)

### Fixed

Many small bugs - refer to commits for details...

### Added

- **Rescaling factor to pre-processed B1**:
To deal with B1 maps that are not in p.u. of the nominal flip angle.
Percent units expected by the hMRI-toolbox, but...
The current BIDS proposal is to have B1 maps scaled so that a value of 1 corresponds to the nominal flip angle.
To be continued...

- **Logging processing messages**:
In order to review and keep track of info, warnings and other messages coming up during data processing, all messages are logged using the hmri_log.m script.
Various options available (pop-up messages, messages logged to the Matlab Command Window, messages saved into a log file).
Improved tracking and readability of the messages, with more explicit descriptions and improved uniformity of the used format.

### Changed

- **Imperfect spoiling correction disabled by default**:
To avoid confusion and resulting mistake (applying correction coefficients to the wrong sequence).
With the publication of "standard MPM protocols", it is likely to have new protocols implemented with sequences other than the customised sequences for which the correction coefficients have been calculated.
Therefore the TR and FA criteria are not sufficient any longer to identify the right set of correction coeficients to be used unambiguously :/...


## [v0.1.2-beta] (released 2018-06-15)

### Fixed

Many small bugs - refer to commits for details...

### Added

- **Single echo VFA approach**:
Implemented with a default minimum number of echoes for R2* calculation of 4.
The number of echoes used is logged and a general warning informs the user
about data interpretation when R2* map derived from only a few echoes.

- **Metadata & processing parameters**:
Bunch of modifications to improve the readability and completeness of the metadata, for each output image as well as for the processing parameters.

- **UNICORT-generated R1 and B1 to calculate PD and MT maps**:
This option has not been thoroughly tested.
Made available under "ADVANCED USERS ONLY": the option to use the R1-UNICORT-derived
B1 map for B1 transmit bias correction in PD and MT maps can only be
enabled by an advanced user implementing customized defaults.

- **Example files**: For defaults customization and toolbox configuration examples.

- **New option to disable the coregistration steps**:
Coregistration shouldn't be disabled but it can be convenient in specific
cases (simulated data, phantom data).
Option made available under "ADVANCED USERS ONLY", i.e. can only be
modified by an advanced user implementing customized defaults.
When enabled (defaults), all input images to hmri_create_MTProt.m (transmit and receive
fields, T1w and MTw images) are coregistered to the PDw average (or TE=0
fit) image (see hmri_create_MTProt.m).

### Changed

- **Options for RF sensitivity bias field correction**:
The available options are now `None`, `Unified Segmentation`(default), `Single` and `Per contrast`.
Beware that the previously implemented `None`option corresponded to the current `Unified Segmentation`!!!
Batch files saved with the previous version are not fully compatible with the current version and must be adapted.


## [v0.1.1-beta] (released 2017-11-14)

### First public beta-version of the hMRI-toolbox.

Modular structure including the following modules:
- `Configure toolbox`
- `DICOM Import`
- `Auto-Reorient`
- `Create hMRI maps`
- `Process hMRI maps`

## [v0.1.0] (released 2016-11-02)

### First stable version of the VBQ toolbox

Considered as mainstream among the various sites taking part in the project
and starting point for the hMRI-toolbox development.
