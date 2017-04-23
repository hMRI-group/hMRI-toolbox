function VO = hmri_pm_segment(InputImage)

% To segment the brain and extract a brain mask in the hMRI Toolbox
% This function replaces pm_segment used in the FieldMap toolbox
% (SPM12/toolbox/FieldMap). Used in hmri_pm_brain_mask.
%
% FORMAT: VO = hmri_pm_segment(P.fname).
%
%========================================================================== 
% Written by Lester Melie-Garcia
% LREN, CHUV, Lausanne, October 2014
%========================================================================== 

[ImagePath,ImageName,ImageExt] = fileparts(InputImage);
Vimg = spm_vol(InputImage);

writeOpts.biascor = 1;
writeOpts.GM  = [0 0 1];
writeOpts.WM  = [0 0 1];
writeOpts.CSF = [0 0 1];
%writeOpts.EXTRA1 = [0 0 1];
%writeOpts.EXTRA2 = [0 0 1];

writeOpts.cleanup = 0;

results   = spm_preproc(Vimg);
[po,pin] = spm_prep2sn(results); %#ok
spm_preproc_write(po,writeOpts);

GMImage  = [ImagePath,filesep,'c1',ImageName,ImageExt];
WMImage  = [ImagePath,filesep,'c2',ImageName,ImageExt];
CSFImage = [ImagePath,filesep,'c3',ImageName,ImageExt];

VO(1).dat = spm_read_vols(spm_vol(GMImage));
VO(2).dat = spm_read_vols(spm_vol(WMImage));
VO(3).dat = spm_read_vols(spm_vol(CSFImage));

delete(GMImage); delete(WMImage); delete(CSFImage);

%% Parameters for spm_preproc.m. We are taking in this case the default parameters
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

%% Parameters coming from old pm_segment
%  flags.fwhm=5;
%  flags.nerode=2;
%  flags.ndilate=4;
%  flags.thresh=0.5;
%  flags.reg = 0.02;
%  flags.graphics=0;
   
end

