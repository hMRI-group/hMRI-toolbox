function [allub1_img] = hmri_create_B1Map_process(ub1_img,ustd_img,vdm_img,fpm_img,b1map_params)
%==========================================================================
% PURPOSE
% To mask, pad and smooth the unwrapped B1 map.
% Part of the hMRI toolbox, B1+ map calculation, EPI SE/STE protocol.
%
% FORMAT [allub1_img] = hmri_create_B1Map_process(ub1_img,ustd_img,vdm_img,fpm_img,b1map_params)
% INPUT ARGUMENTS (all filenames are 1-element cell arrays)
% - ub1_img         filename of the unwrapped B1 map  
% - ustd_img        filename of the unwrapped SD (error) image 
% - vdm_img         filename of the matched VDM image 
% - fpm_img         filename of the phase-unwrapped regularised field map (Hz)
% - b1map_params    effective parameters defined for the current B1 map
%                       calculation, including B1 processing parameters.
% OUTPUT ARGUMENTS
% - allub1_img{1}   Masked padded unwarped B1 map
% - allub1_img{2}   Smoothed masked padded unwarped B1 map
%==========================================================================
% Written by Antoine Lutti, based on Chloe Hutton's code and the FieldMap
% toolbox. 

if nargin ~= 5
  error('Enter unwarped B1, SD(error) map, matched VDM file, Hz fieldmap and B1map parameters');
end

%--------------------------------------------------------------------------
% retrieves default parameters for B1 processing:
%--------------------------------------------------------------------------
b1proc = b1map_params.b1proc;

%--------------------------------------------------------------------------
% Make a thresholded mask of the B1 STD map
%--------------------------------------------------------------------------
   
stdmask = spm_read_vols(spm_vol(ustd_img{1}(1,:)))<=b1proc.SDTHRESH;
   
%--------------------------------------------------------------------------
% Make a thresholded mask of the Hz field 
%--------------------------------------------------------------------------
Vfpm = spm_vol(fpm_img{1}(1,:));
hzmask = spm_read_vols(Vfpm);
hzmask = hzmask<b1proc.HZTHRESH & hzmask>-b1proc.HZTHRESH;

% Erode the Hz mask a bit
nerode = b1proc.ERODEB1;
Vb1 = spm_vol(ub1_img{1}(1,:));
pxs = sqrt(sum(Vb1.mat(1:3,1:3).^2)); % Voxel resolution
npxs = pxs/pxs(1);
ks = floor((nerode+1)*[1 1 1]./npxs);
if ~mod(ks(1),2); ks(1) = ks(1)+1;end
if ~mod(ks(2),2); ks(2) = ks(2)+1;end
if ~mod(ks(3),2); ks(3) = ks(3)+1;end
kernel = ones(ks);
if size(size(kernel),2)~=3
       kk(:,:,1)=kernel;
       kk(:,:,2)=kernel;
       kernel=kk;
       clear kk;
end
hzmask = spm_erode(double(hzmask),kernel);

% Resample mask in space of b1map:
[x,y,z] = ndgrid(1:Vb1.dim(1),1:Vb1.dim(2),1:Vb1.dim(3));

% Need V.mat from the matched VDM file:
Vvdm = spm_vol(vdm_img{1}(1,:));
tM = inv(Vb1.mat\Vvdm.mat);

x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
xyz2 = [x2(:) y2(:) z2(:)];
rhzmask = double(reshape(spm_sample_vol(hzmask,xyz2(:,1),xyz2(:,2),xyz2(:,3),1),Vb1.dim(1:3))>0);

% Multiply rhz mask by the the std mask
fullmask = rhzmask.*stdmask; 
nerode = b1proc.PADB1;
ks = floor((nerode+1)*[1 1 1]./npxs);
if ~mod(ks(1),2); ks(1) = ks(1)+1;end
if ~mod(ks(2),2); ks(2) = ks(2)+1;end
if ~mod(ks(3),2); ks(3) = ks(3)+1;end
kernel = ones(ks);
if size(size(kernel),2)~=3
       kk(:,:,1) = kernel;
       kk(:,:,2) = kernel;
       kernel = kk;
       clear kk;
end

%--------------------------------------------------------------------------
% Pad the masked B1 map
%--------------------------------------------------------------------------
b1map = spm_read_vols(Vb1);
[padb1map,padmask] = pm_pad(b1map.*fullmask,fullmask,kernel); %#ok<NASGU>
Vpadb1 = Vb1;
Vpadb1.pinfo(1:2) = Inf;
[pth,fname,ext] = fileparts(Vb1.fname);
Vpadb1.fname = fullfile(pth,['m' fname ext]);
Vpadb1.descrip = 'Masked padded unwarped B1 map';
allub1_img{1} = spm_write_vol(Vpadb1,padb1map);

%--------------------------------------------------------------------------
% Smooth padded B1 map
%--------------------------------------------------------------------------
spadb1map = zeros(size(padb1map));
smth = b1proc.B1FWHM./pxs;
spm_smooth(padb1map,spadb1map,smth);
Vspadb1 = Vb1;
Vspadb1 = rmfield(Vspadb1,'pinfo');
[pth,fname,ext] = fileparts(Vb1.fname);
Vspadb1.fname = fullfile(pth,['sm' fname ext]);
Vspadb1.descrip = sprintf('Smoothed (%dmm) masked padded unwarped B1 map',b1proc.B1FWHM(1));
allub1_img{2} = spm_write_vol(Vspadb1,spadb1map);

end