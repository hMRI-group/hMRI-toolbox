function vbq_B1map_v2long(P,Q,T1)
% Calculation of B1 mapping data 3D EPI spin echo (SE) and stimulated (STE) echo images (see Jiru and Klose MRM 2006).
% Corresponding scanning protocol/sequence: al_B1mapping_v2a and
% al_B1mapping_v2b
% Input: 15 pairs of (SE, STE) images for B1 map calculation and 3 images
% for B0 map calculation.
% This macro calls the functions B1Map_unwarp and B1Map_process for
% correction of image distortions, padding and smoothing of the images.
% Output:
%     - distorted B1 (B1map_----) and error (SDmap_----) maps
%     - undistorted B1 (uB1map_----) and error (uSDmap_----) maps
%     - undistorted, masked and padded B1 maps (muB1map_---------)
%     - undistorted, masked, padded and smoothed B1 maps (smuB1map_---------) i.e. FULLY PROCESSED
% At each voxel, this macro selects the 5 pairs of (SE,STE image) (out of
% 15) with maximum signal amplitude in the SE images.
% The sum of square image of all SE images is created (SumOfSq) and
% undistorted (uSumOfSq) for coregistration of the B1 map to an anatomical dataset

% $Id$

% if nargin == 3: B1 mapping for QA: T1 is given
% else:
if nargin < 3 
    T1 = 1192; %ms
    if nargin < 2
        P = spm_select(Inf,'image','Select images for B1 maps');
        Q = spm_select(Inf,'image','Select images for B0 map');
    end
end

beta=[135:-5:65];
TM = 33.24;
eps=0.0001;
Nonominalvalues=5;

disp('----- Calculation of B1 map -----');

V=spm_vol(P);
n   = numel(V);
AllPossB1   = zeros([V(1).dim(1:2) n]);%All Possible values of the local B1 value calculated from alocal=acos() and 180-alocal
OptB1 = zeros([V(1).dim(1:2) 2*Nonominalvalues]);%Possible B1 estimates (alocal and 180-alocal) with maximum SE signal 
B1local = zeros([V(1).dim(1:3)]);
SDlocal = zeros([V(1).dim(1:3)]);
Ssq_matrix=sqrt(sum(spm_read_vols(V(1:2:end-1)).^2,4));

%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',V(1).dim(3),'B1 map fit','planes completed');

%-Loop over planes computing result Y
%-----------------------------------------------------------------------

corr_fact = exp(TM/T1);
for p = 1:V(1).dim(3),%loop over the partition dimension of the data set
    SEsignal=zeros(V(1).dim(1), V(1).dim(2), n/2);
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
    M = inv(B*inv(V(1).mat)*V(1).mat);
    for i = 1:n/2
        AllPossB1(:,:,i)  = real(acos(corr_fact*spm_slice_vol(V((i-1)*2+2),M,V(1).dim(1:2),0)./(spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0)+eps))/pi*180/beta(i)); % nearest neighbor interpolation
        AllPossB1(:,:,n/2+i)  = 180/beta(i) - AllPossB1(:,:,i);
        SEsignal(:,:,i)=spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0);
    end
    [~, indexes]=sort(SEsignal,3);%The optimal B1 estimates are obtained when the signal amplitude in the SE image is maximum
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
            OptB1(x_nr,y_nr,1:Nonominalvalues)=AllPossB1(x_nr,y_nr,indexes(x_nr,y_nr,n/2-(Nonominalvalues-1):end));
            OptB1(x_nr,y_nr,Nonominalvalues+1:2*Nonominalvalues)=AllPossB1(x_nr,y_nr,indexes(x_nr,y_nr,n/2-(Nonominalvalues-1):end)+n/2);
        end
    end

    OptB1=sort(real(OptB1), 3); % take the real value due to noise problems
    Y_sd   = zeros([V(1).dim(1:2) (Nonominalvalues+1)]);
    Y_mn   = zeros([V(1).dim(1:2) (Nonominalvalues+1)]);
    for i = 1:(Nonominalvalues+1)
       Y_sd(:,:,i) = std(OptB1(:,:,i:(i+Nonominalvalues-1)), [], 3);
       Y_mn(:,:,i) = mean(OptB1(:,:,i:(i+Nonominalvalues-1)), 3);
    end
    [~, min_index] = min(Y_sd,[],3);%!!min_index is a 2D array. Size given by resolution along read and phase directions
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
             B1local(x_nr,y_nr,p) = Y_mn(x_nr,y_nr, min_index(x_nr,y_nr));%Y_ab is the relative flip angle value averaged over the n flip angles (determined by minizing the SD i.e. keeping the most uniform relative flip angle values)
             SDlocal(x_nr,y_nr,p) = Y_sd(x_nr,y_nr, min_index(x_nr,y_nr));%Y_cd is the corresponding standard deviation between the relative flip angle values
        end
    end
    spm_progress_bar('Set',p);
end

V_save=struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','B1 map [%]');
V_save.fname=fullfile(spm_str_manip(V_save.fname,'h'),['B1map_' spm_str_manip(V_save.fname,'t')]);
V_save = spm_write_vol(V_save,B1local*100);

W_save=struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SD [%]');
W_save.fname=fullfile(spm_str_manip(W_save.fname,'h'),['SDmap_' spm_str_manip(W_save.fname,'t')]);
W_save = spm_write_vol(W_save,SDlocal*100);

X_save=struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SE SSQ matrix');
X_save.fname=fullfile(spm_str_manip(X_save.fname,'h'),['SumOfSq.' spm_str_manip(X_save.fname,'e')]);
X_save = spm_write_vol(X_save,Ssq_matrix);
 
pm_defaults;
pm_defs=pm_def;
pm_defs.SHORT_ECHO_TIME=10;
pm_defs.LONG_ECHO_TIME=12.46;
pm_defs.blipdir=1;
pm_defs.TOTAL_EPI_READOUT_TIME=540e-3*24;
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

anat_img = {char(X_save.fname)};
other_img{1}=char(V_save.fname);
other_img{2}=char(W_save.fname);

[fmap_img,unwarp_img] = B1Map_unwarp(fm_imgs,anat_img,other_img,pm_defs);
uanat_img{1}=unwarp_img{1}.fname;
ub1_img{1}=unwarp_img{2}.fname;
ustd_img{1}=unwarp_img{3}.fname;
fpm_img{1}=fmap_img{1};
vdm_img{1}=fmap_img{2};
[allub1_img]=B1Map_process(uanat_img,ub1_img,ustd_img,vdm_img,fpm_img,pm_defs);
