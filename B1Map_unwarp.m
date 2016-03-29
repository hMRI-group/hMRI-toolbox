function [fmap_img,unwarp_img]=B1Map_unwarp(fm_imgs,anat_img,other_img,pm_defs)
% Ensure that your version of Fieldmap in 'Fieldmapxy' is the right
% version.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id$

if nargin < 4
  error('Enter field map images, distorted anatomical image, images to unwarp and defaults');
end

IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% Define parameters for fieldmap creation

if ~isfield(pm_defs,'et')
   IP.et{1} = pm_defs.SHORT_ECHO_TIME;
   IP.et{2} = pm_defs.LONG_ECHO_TIME;
else
   IP.et{1} = pm_defs.et{1};
   IP.et{2} = pm_defs.et{2};
end

if ~isfield(pm_defs,'maskbrain')
   IP.maskbrain = pm_defs.MASKBRAIN;
else
   IP.maskbrain = pm_defs.maskbrain;
end
     
% Set parameters for unwrapping
if ~isfield(pm_defs,'uflags')
   IP.uflags.iformat = pm_defs.INPUT_DATA_FORMAT;
   IP.uflags.method = pm_defs.UNWRAPPING_METHOD;
   IP.uflags.fwhm = pm_defs.FWHM;
   IP.uflags.pad = pm_defs.PAD;
   IP.uflags.ws = pm_defs.WS;
   IP.uflags.etd = pm_defs.LONG_ECHO_TIME - pm_defs.SHORT_ECHO_TIME;     
else
   IP.uflags.iformat = pm_defs.uflags.iformat;
   IP.uflags.method = pm_defs.uflags.method;
   IP.uflags.fwhm = pm_defs.uflags.fwhm;
   IP.uflags.pad = pm_defs.uflags.pad;
   IP.uflags.ws = pm_defs.uflags.ws;
   IP.uflags.etd = pm_defs.uflags.etd;      
end

% Set parameters for brain extraction
if ~isfield(pm_defs,'mflags')
   IP.mflags.template=pm_defs.MFLAGS.TEMPLATE;
   IP.mflags.fwhm=pm_defs.MFLAGS.FWHM; 
   IP.mflags.nerode=pm_defs.MFLAGS.NERODE;
   IP.mflags.ndilate=pm_defs.MFLAGS.NDILATE;
   IP.mflags.thresh=pm_defs.MFLAGS.THRESH;
else
   IP.mflags.template=pm_defs.mflags.template;
   IP.mflags.fwhm=pm_defs.mflags.fwhm; 
   IP.mflags.nerode=pm_defs.mflags.nerode;
   IP.mflags.ndilate=pm_defs.mflags.ndilate;
   IP.mflags.thresh=pm_defs.mflags.thresh;
end

% Get FieldMap parameters 
if ~isfield(pm_defs,'ajm')
   IP.ajm = pm_defs.DO_JACOBIAN_MODULATION;
else
   IP.ajm = pm_defs.ajm;
end
if ~isfield(pm_defs,'blipdir')
   IP.blipdir = pm_defs.K_SPACE_TRAVERSAL_BLIP_DIR;
else
   IP.blipdir = pm_defs.blipdir;
end
if ~isfield(pm_defs,'tert')
   IP.tert = pm_defs.TOTAL_EPI_READOUT_TIME;
else
   IP.tert = pm_defs.tert;
end 
if ~isfield(pm_defs,'epifm')
   IP.epifm = pm_defs.EPI_BASED_FIELDMAPS;
else
   IP.epifm = pm_defs.epifm; 
end

% Clear any old handles etc
IP.fm = [];
IP.vdm = [];
IP.jim = [];
IP.pP = [];
IP.epiP = [];
IP.uepiP = [];
IP.vdmP = [];
ID = cell(4,1);

%----------------------------------------------------------------------
% Load measured field map data - phase and magnitude or real and imaginary
%----------------------------------------------------------------------

if ischar(fm_imgs)
   fm_imgs = spm_vol(fm_imgs); 
end
n_fms = length(fm_imgs);
switch n_fms
case 4  % real, imaginary pairs
   for i = 1:n_fms
      IP.P{i} = spm_vol(fm_imgs(i)); 
   end
case 2  % precalculated phase map and magnitude image
   IP.P{1} = spm_vol(fm_imgs(1));
   IP.P{2} = spm_vol(fm_imgs(2));
case 1  % precalculated and unwrapped Hz map    
   IP.pP = spm_vol(fm_imgs);
   IP.fm.fpm = spm_read_vols(IP.pP);
   IP.fm.jac = pm_diff(IP.fm.fpm,2);
   if isfield(IP,'P') & ~isempty(IP.P{1})
      IP.P = cell(1,4); 
   end
   if isfield(pm_defs,'magfieldmap')
      IP.fmagP=pm_defs.magfieldmap;
   end
otherwise 
   error('Funny number of input fieldmap images')
end

if ~isempty(IP.P{1})
%----------------------------------------------------------------------
% Create field map (in Hz) - this routine calls the unwrapping
%----------------------------------------------------------------------
   IP.fm = FieldMap('CreateFieldMap',IP);
%----------------------------------------------------------------------
% Write out field map
% Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------  
   fmap_img{1}=FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Smoothed phase map');
end
%----------------------------------------------------------------------
% Convert Hz to voxels and write voxel displacement map 
% Outputs -> vdm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

[IP.vdm, IP.vdmP]=FieldMap('FM2VDM',IP);

%----------------------------------------------------------------------
% Match voxel displacement map to distorted anatomical image
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------
% To DO!!!

%----------------------------------------------------------------------
% Match voxel displacement map to distorted anatomical image and 
% then unwarp the distorted anatomical image
%----------------------------------------------------------------------
IP.epiP = spm_vol(anat_img{1}(1,:));
if isfield(pm_defs, 'match_vdm')
   if pm_defs.match_vdm
          IP.vdmP = FieldMap('MatchVDMxy',IP);
   end
end
IP.uepiP = FieldMap('UnwarpEPIxy',IP);
unwarp_info=sprintf('Unwarped image:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
unwarp_img{1}=FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);
%----------------------------------------------------------------------
% Unwarp the other images, it is assumed they are in same space as
% anatomical
%----------------------------------------------------------------------
for filenum=1:numel(other_img)
    
   IP.epiP = spm_vol(other_img{filenum}(1,:));
   IP.uepiP = FieldMap('UnwarpEPIxy',IP);
    
   %----------------------------------------------------------------------
   % Write unwarped EPI 
   % Outputs -> uNAME-OF-EPI.img
   %----------------------------------------------------------------------
   unwarp_info=sprintf('Unwarped image:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
   unwarp_img{1+filenum}=FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);
   fmap_img{2}=IP.vdmP;
   %IPcell{1}=IP;

end








