function [fmap_img,unwarp_img] = hmri_create_B1Map_unwarp(fmfnam, anatfnam, otherfnam, b1map_params)
%==========================================================================
% PURPOSE
% For B0 undistortion of EPI-based B1 maps, part of the hMRI toolbox.
%
% FORMAT [fmap_img,unwarp_img] = hmri_create_B1Map_unwarp(magfnam, phasefnam, anatfnam, otherfnam, b1map_params)
% INPUT ARGUMENTS
% - fmfnam          filename of the fieldmap images (must be char) 
% - anatfnam        filename of the image to be undistorted ("anatomical" reference, i.e. SSQ image)   
% - otherfnam       filenames of the other images to undistort (B1 and SD maps)
% - b1map_params    effective parameters defined for the current B1 map
%                       calculation, including B0 unwarping parameters.
% OUTPUT ARGUMENTS
% - fmap_img{1} = IP.fpm = phase-unwrapped regularised field map (Hz)
% - fmap_img{2} = IP.vdmP = vdm5_* file name
% - unwarp_img{1} = IP.uepiP = unwrapped anat image
% - unwarp_img{2,3,...} = unwrapped other image{1,2,...}
%==========================================================================
% Written by Antoine Lutti, based on Chloe Hutton's code and the FieldMap
% toolbox. 

if nargin < 4
  error('Enter field map images, distorted anatomical image, images to unwarp and defaults');
end

% initialise parameters (default values)
pm_defaults;
pm_defs = pm_def;
IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% fill in with effective parameters for fieldmap undistortion
IP.et{1} = b1map_params.b0acq.shortTE;
IP.et{2} = b1map_params.b0acq.longTE;
IP.maskbrain = b1map_params.b1proc.b0maskbrain;
IP.ajm = pm_defs.DO_JACOBIAN_MODULATION;
IP.blipdir = b1map_params.b1acq.blipDIR;
IP.tert = b1map_params.b1acq.tert;
%IP.epifm = 0; % already in the default IP
%IP.pedir = 2; % already in the default IP ([1] A>>P, [2] R>>L) 

% default unwarp parameters (uflags) and a few derived from the metadata 
IP.uflags.iformat = b1map_params.b0acq.iformat; % input format = 'PM'
IP.uflags.method = pm_defs.UNWRAPPING_METHOD;
IP.uflags.fwhm = pm_defs.FWHM;
IP.uflags.pad = pm_defs.PAD;
IP.uflags.ws = pm_defs.WS;
IP.uflags.etd = IP.et{2}-IP.et{1};

% default brain extraction parameters (mflags) 
IP.mflags.template = pm_defs.MFLAGS.TEMPLATE;
IP.mflags.fwhm = pm_defs.MFLAGS.FWHM; 
IP.mflags.nerode = pm_defs.MFLAGS.NERODE;
IP.mflags.ndilate = pm_defs.MFLAGS.NDILATE;
IP.mflags.thresh = pm_defs.MFLAGS.THRESH;

% Clear any old handles etc
IP.fm = [];
IP.vdm = [];
IP.jim = [];
IP.pP = [];
IP.epiP = [];
IP.uepiP = [];
IP.vdmP = [];

%--------------------------------------------------------------------------
% Load measured field map data - phase and magnitude or real and imaginary
% NOTE: only the magnitude + presubtracted phase option has been used and
% tested at this stage, therefore the other options are commented out.
%--------------------------------------------------------------------------
n_fms = size(fmfnam,1);
switch n_fms
%     case 4  % real, imaginary pairs
%         for i = 1:n_fms
%             IP.P{i} = spm_vol(fmfnam(i,:));
%         end
    case 2  % precalculated phase map and magnitude image
        % first rescale phase map between +/-pi
        scphase = FieldMap('Scale', fmfnam(1,:));
        try
            movefile(scphase.fname, fullfile(b1map_params.outpath, spm_file(scphase.fname,'filename')));
            scphase.fname = fullfile(b1map_params.outpath, spm_file(scphase.fname,'filename'));
        catch % will fail if trying to movefile onto itself
        end
        IP.P{1} = spm_vol(scphase.fname);
        IP.P{2} = spm_vol(fmfnam(2,:));
%     case 1  % precalculated and unwrapped Hz map
%         IP.pP = spm_vol(fmfnam);
%         IP.fm.fpm = spm_read_vols(IP.pP);
%         IP.fm.jac = pm_diff(IP.fm.fpm,2);
%         if isfield(IP,'P') & ~isempty(IP.P{1})
%             IP.P = cell(1,4);
%         end
%         if isfield(pm_defs,'magfieldmap')
%             IP.fmagP=pm_defs.magfieldmap;
%         end
    otherwise
        error('Funny number of input fieldmap images')
end

if ~isempty(IP.P{1})
    %----------------------------------------------------------------------
    % Create field map (in Hz) - this routine calls the unwrapping
    %----------------------------------------------------------------------
    IP.fm = hmri_create_FieldMap('CreateFieldMap',IP);
    
    % TL: move created maps to outpath
    if IP.maskbrain==1
        % move bmask file (brain mask) and m file (magnitude image normalised (?))
        try
            bmaskfnam = spm_file(IP.P{2}.fname, 'prefix','bmask');
            movefile(bmaskfnam, fullfile(b1map_params.outpath, spm_file(bmaskfnam,'filename')));
            mfnam = spm_file(IP.P{2}.fname, 'prefix','m');
            movefile(mfnam, fullfile(b1map_params.outpath, spm_file(mfnam,'filename')));
        catch % will fail if trying to movefile onto itself
        end
    end
    %----------------------------------------------------------------------
    % Write out IP.fm.fpm = phase-unwrapped, regularised field map (Hz)
    % Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.nii
    % NB: no need to movefile since it is written in same directory as
    % scphase which had been moved already.
    %----------------------------------------------------------------------
    fmap_img{1} = FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Smoothed phase map');
end

%--------------------------------------------------------------------------
% Convert Hz to voxels and write voxel displacement map 
% Outputs -> vdm5_NAME-OF-FIRST-INPUT-IMAGE.nii
% NB: no need to movefile since it is written in same directory as
% scphase which had been moved already.
%--------------------------------------------------------------------------
[IP.vdm, IP.vdmP] = FieldMap('FM2VDM',IP);

%--------------------------------------------------------------------------
% Match voxel displacement map to distorted anatomical image
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.nii
%--------------------------------------------------------------------------
% To DO!!!

%--------------------------------------------------------------------------
% Match voxel displacement map to distorted anatomical image and 
% then unwarp the distorted anatomical image
%--------------------------------------------------------------------------
IP.epiP = spm_vol(anatfnam);
if b1map_params.b1proc.match_vdm
    IP.vdmP = FieldMap('MatchVDMxy',IP); 
    % write forward warped magnitude image (i.e. magnitude image from
    % fieldmap acquisition, coregistered & resliced to anat image)
    % wfmag_NAME-OF-ANAT-IMAGE.nii. Result is written directly in output
    % directory, no need to movefile. Also coregister (not reslice) vdm5_*
    % image with it (IP.vdmP.fname = vdm5_* file name). 
end
IP.uepiP = FieldMap('UnwarpEPIxy',IP);
unwarp_info = sprintf('Unwarped image:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
% write uNAME-OF-ANAT-IMAGE.nii in output directory, no need to movefile.
unwarp_img{1} = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);

% set second fmap_img output to the coregistered VDM
fmap_img{2} = IP.vdmP;


%----------------------------------------------------------------------
% Unwarp the other images, assumed to be in same space as anatomical
%----------------------------------------------------------------------
for filenum = 1:numel(otherfnam)
    
   IP.epiP = spm_vol(otherfnam{filenum}(1,:));
   IP.uepiP = FieldMap('UnwarpEPIxy',IP);
    
   %----------------------------------------------------------------------
   % Write unwarped EPI 
   % Outputs -> uNAME-OF-EPI.nii
   %----------------------------------------------------------------------
   unwarp_info = sprintf('Unwarped image:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
   OrigMat = spm_read_vols(IP.epiP);
   NewMat = IP.uepiP.dat;
   % AL - 05/10/2016
   % In the unwarped images, replace voxels set to zero by the Fieldmap
   % toolbox by their original values in the warped data. It's a fudge but:
   % 1. Not expected to do any harm in normal circumstances (voxels set to
   %    0 are outside the brain) 
   % 2. Helps when fieldmap data doesn't quite cover the whole brain due to
   %    FOV restrictions: in this case, the voxels in the B1 (and SD) maps
   %    outside the fieldmap FOV are set to 0
   NewMat(NewMat==0) = OrigMat(NewMat==0);
   IP.uepiP.dat = NewMat;
   unwarp_img{1+filenum} = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info); %#ok<AGROW>

end

end