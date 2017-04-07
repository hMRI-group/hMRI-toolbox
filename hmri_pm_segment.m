function VO = hmri_pm_segment(InputImage)

% This function was developed to substitute the function pm_segment used in
% FieldMap toolbox : ...\spm12b\toolbox\FieldMap\pm_segment.m
% 
%This functions replace line:
% VO=pm_segment(P.fname,flags.template,seg_flags) in line 48 of function 'pm_brain_mask.m'
%
% It will be:  VO = pm_segment_lren(P.fname). It could be improved
% customizing the segmentation parameters. In this case the default values are taken.
%
%% Lester Melie-Garcia
% LREN, CHUV. 
% Lausanne, October 2nd, 2014

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

