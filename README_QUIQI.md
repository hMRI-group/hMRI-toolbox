
# Inserting an index of motion degradation into the statistical analysis of MRI data (QUIQI): Toolbox Implementation

### Authors

- Antoine Lutti, 2021
  Laboratory for Neuroimaging Research
  Lausanne University Hospital & University of Lausanne, Lausanne, Switzerland

- Nadège Corbin, 2021
  Centre de Résonance Magnétique des Systèmes Biologiques, Bordeaux, France

Copyright (C) 2021 Laboratory for Neuroimaging Research

## REFERENCE
Lutti et al. ‘Restoring statistical validity in group analyses of motion-corrupted MRI data’, Human Brain Mapping, 2022. https://doi.org/10.1002/hbm.25767

## INTRODUCTION
QUIQI is a method that accounts for the degradation of data quality due to motion in the analysis of MRI data. This method is described in the above-mentioned scientific article. QUIQI was originally implemented using custom-made matlab scripts and functionalities provided by the SPM software (https://www.fil.ion.ucl.ac.uk/spm/). These scripts are available here: (https://github.com/LREN-physics/hMRI-toolbox/tree/quiqi_v1.0). A subset of the data analysed in the original article can be found here: doi:10.5281/zenodo.4647081.
In order to help users implement this method for their own analysis, we provide here a GUI-based implementation of QUIQI built into the hMRI toolbox, dedicated to the analysis of neuroimaging data (https://hmri-group.github.io/hMRI-toolbox/). At the time of writing, the version of the hMRI toolbox with the QUIQI implementation is separate from the main toolbox branch. We expect merging of this version into the main branch to take place in the near future, following further development.

## REQUIREMENTS
A running version of Matlab (https://ch.mathworks.com/products/matlab.html) and SPM (https://www.fil.ion.ucl.ac.uk/spm/) are pre-requisites for running these analyses. Please refer to https://hmri-group.github.io/hMRI-toolbox/ for help in installing and using the hMRI toolbox.

## CONTENT
The implementation of QUIQI within the hMRI toolbox is based on two modules: ‘QUIQI BUILD’ and ‘QUIQI CHECK’.

### 1. QUIQI BUILD
The ‘QUIQI BUILD’ module lies between the usual ‘Factorial design’ and ‘Model estimation’ steps of an SPM image analysis (see screenshot below). The corresponding fields are:
- SPM.mat: SPM matrix associated with the analysis
- MDI powers: powers of the Motion Degradation Index used to compute the basis functions for the ReML estimation of the noise covariance matrix.
- MDI type: type of the object that contains the values of the Motion Degradation Index. If ‘MDI matrix’ is selected, the MDI values may be pasted into the GUI from e.g. an Matlab variable or a table. The MDI values can also be read-in directly from a json file (‘MDIjsonfile).

![QUIQI_build](figs\QUIQI_build.png)

### 2. QUIQI CHECK
QUIQI CHECK lies at the very end of image analysis (see below, left). While this module is not necessary to correct for motion degradation effects in the MRI data, it is useful to assess the degree of remaining motion degradation after correction. On the model of Lutti et al., this module plots a histogram of the distribution of the image residual samples against the results from the fitting of these residuals as a polynomial function of the MDI. A high R2 value highlights heteroscedasticity of the noise distribution (see example below, with and without the use of QUIQI). Note that when multiple sets of MDI values are used for QUIQI (i.e. for the analysis of MRI data computed across multiple raw image volumes) all sets of MDI values are used for the polynomial fitting. Also, the maximum power of the polynomial fit can be set from the user interface (field ‘Fit power’).

<img src="figs\QUIQI_check1.png" alt="QUIQI_check1" style="zoom:80%;" />
<img src="figs\QUIQI_check2.png" alt="QUIQI_check2" style="zoom:67%;" /><img src="figs\QUIQI_check3.png" alt="QUIQI_check3" style="zoom: 67%;" />