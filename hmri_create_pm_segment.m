function VO = hmri_create_pm_segment(InputImage)
% To segment the brain and extract a brain mask in the hMRI Toolbox
% This function replaces pm_segment used in the FieldMap toolbox
% (SPM12/toolbox/FieldMap). Used in hmri_create_pm_brain_mask.
%
% FORMAT: VO = hmri_create_pm_segment(P.fname).
%
%========================================================================== 
% Written by Lester Melie-Garcia
% LREN, CHUV, Lausanne, October 2014
% Adapted by Evelyne Balteau, CRC, Liège, September 2017...
% to use unified segmentation for "forward" compatibility!
%========================================================================== 

% use unified segmentation instead of OldSegment
% use default parameters:
job = hmri_get_defaults('segment');
job.channel.vols = {InputImage};
% except for:
for ctis=4:length(job.tissue)
    job.tissue(ctis).native = [0 0]; % no need to write c4, c5, ...
end
job.channel.write = [0 0]; % no need to write BiasField nor BiasCorrected volume
% run segmentation
segm_output = spm_preproc_run(job);

GMImage  = segm_output.tiss(1).c{1};
WMImage  = segm_output.tiss(2).c{1};
CSFImage = segm_output.tiss(3).c{1};

VO(1).dat = spm_read_vols(spm_vol(GMImage));
VO(2).dat = spm_read_vols(spm_vol(WMImage));
VO(3).dat = spm_read_vols(spm_vol(CSFImage));

delete(GMImage); delete(WMImage); delete(CSFImage);

% % Parameters for spm_preproc.m. We are taking in this case the default parameters
%  opts - options
%  opts.tpm      - n tissue probability images for each class
%  opts.ngaus    - number of Gaussians per class (n+1 classes)
%  opts.warpreg  - warping regularisation
%  opts.warpco   - cutoff distance for DCT basis functions
%  opts.biasreg  - regularisation for bias correction
%  opts.biasfwhm - FWHM of Gausian form for bias regularisation
%  opts.regtype  - regularisation for affine part
%  opts.fudge    - a fudge factor
%  opts.msk      - unused

% % Parameters coming from old pm_segment
%  flags.fwhm=5;
%  flags.nerode=2;
%  flags.ndilate=4;
%  flags.thresh=0.5;
%  flags.reg = 0.02;
%  flags.graphics=0;
   
end

