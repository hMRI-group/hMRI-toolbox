function B1map_rect700(P, Q) % MFC: added arguments P & Q to mirror B1map_v2
% Calculation of B1 mapping data 3D EPI spin echo (SE) and stimulated (STE) echo images (see Jiru and Klose MRM 2006).
% Corresponding scanning protocol/sequence: al_B1mapping_rect700
% Inputs: SE and STE images for B1 map calculation and the 3 data images
% for B0 map calculation.
% This macro calls the functions B1Map_unwarp and B1Map_process for
% correction of image distortions, padding and smoothing of the images.
% Output:
%     - distorted B1 (B1map_----) and error (SDmap_----) maps
%     - undistorted B1 (uB1map_----) and error (uSDmap_----) maps
%     - undistorted, masked and padded B1 maps (muB1map_---------)
%     - undistorted, masked, padded and smoothed B1 maps (smuB1map_---------) i.e. FULLY PROCESSED
% IN THIS VERSION, THE B1 MAPS ARE CALCULATED BASED ON ALL 5 PAIRS OF SE/STE IMAGES

% $Id$

beta=[80:5:100];
TM = 33.53;%
T1 = 1192; %ms
eps=0.0001;

if nargin < 2
%     P = spm_select(Inf,'-0..nii','Select the SE/STE images'); 
    P = spm_select(Inf,'image','Select images for B1 maps'); 
    Q = spm_select(Inf,'image','Select images for B0 map'); 
end
disp('----- Calculation of B1 map -----');

V=spm_vol(P);
n   = numel(V);
Y_tmp   = zeros([V(1).dim(1:2) n]);
Y_ab = zeros([V(1).dim(1:3)]);
Y_cd = zeros([V(1).dim(1:3)]);
%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',V(1).dim(3),'B1 map fit','planes completed');

%-Loop over planes computing result Y
%-----------------------------------------------------------------------

corr_fact = exp(TM/T1);
for p = 1:V(1).dim(3),%loop over the partition dimension of the data set
	B = spm_matrix([0 0 -p 0 0 0 1 1 1]);

	for i = 1:n/2
		M = inv(B*inv(V(1).mat)*V(1).mat);
        Y_tmp(:,:,((i-1)*2+1))  = real(acos(corr_fact*spm_slice_vol(V((i-1)*2+2),M,V(1).dim(1:2),0)./(spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0)+eps))/pi*180/beta(i)); % nearest neighbor interpolation
        Y_tmp(:,:,((i-1)*2+2))  = 180/beta(i) - Y_tmp(:,:,((i-1)*2+1));
    end
   
    Y_tmp = sort(real(Y_tmp), 3); % take the real value due to noise problems
    
    Y_sd   = zeros([V(1).dim(1:2) (n/2+1)]);
    Y_mn   = zeros([V(1).dim(1:2) (n/2+1)]);
    for i = 1:(n/2+1)
       Y_sd(:,:,i) = std(Y_tmp(:,:,i:(i+n/2-1)), [], 3);
       Y_mn(:,:,i) = mean(Y_tmp(:,:,i:(i+n/2-1)), 3);
    end
    [dummy, min_index] = min(Y_sd,[],3);%!!min_index is a 2D array. Size given by resolution along read and phase directions
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
             Y_ab(x_nr,y_nr,p) = Y_mn(x_nr,y_nr, min_index(x_nr,y_nr));%Y_ab is the relative flip angle value averaged over the n flip angles (determined by minizing the SD i.e. keeping the most uniform relative flip angle values)
             Y_cd(x_nr,y_nr,p) = Y_sd(x_nr,y_nr, min_index(x_nr,y_nr));%Y_cd is the corresponding standard deviation between the relative flip angle values
        end
    end
    spm_progress_bar('Set',p);
end

V_save=struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','B1 map [%]');
[p,name,e] = fileparts(V_save.fname);
P_B1 = fullfile(p,['B1map_' name e]);
V_save.fname = P_B1;
V_save = spm_write_vol(V_save,Y_ab*100);

W_save=struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SD [%]');
[p,name,e] = fileparts(W_save.fname);
P_SDB1 = fullfile(p,['SDmap_' name e]);
W_save.fname = P_SDB1;
W_save = spm_write_vol(W_save,Y_cd*100);
 
pm_defaults;
pm_defs=pm_def;
pm_defs.SHORT_ECHO_TIME=10;
pm_defs.LONG_ECHO_TIME=12.46;
pm_defs.blipdir=1;
pm_defs.TOTAL_EPI_READOUT_TIME=500e-3*48;
pm_defs.MASKBRAIN=1;
pm_defs.STDTHRESH=5;
pm_defs.HZTHRESH=110;
pm_defs.ERODEB1=1;
pm_defs.PADB1=3 ;
pm_defs.B1FWHM=8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
pm_defs.match_vdm=1;

mag1=spm_vol(Q(1,:));
mag=mag1.fname;
phase1=spm_vol(Q(3,:));
phase=phase1.fname;
scphase=FieldMap('Scale',phase);
fm_imgs=str2mat(scphase.fname,mag);

anat_img1=spm_vol(P(5,:));
[path,name,e] = fileparts(anat_img1.fname);
anat_img = {strcat(path,filesep,name,e)};
other_img{1}=char(V_save.fname);
other_img{2}=char(W_save.fname);

[fmap_img,unwarp_img] = B1Map_unwarp(fm_imgs,anat_img,other_img,pm_defs);
uanat_img{1}=unwarp_img{1}.fname;
ub1_img{1}=unwarp_img{2}.fname;
ustd_img{1}=unwarp_img{3}.fname;
fpm_img{1}=fmap_img{1};
vdm_img{1}=fmap_img{2};
[allub1_img]=B1Map_process(uanat_img,ub1_img,ustd_img,vdm_img,fpm_img,pm_defs);
