function P_trans = run_vbq_b1map(jobsubj)

%% Processing of the B1 maps for B1 bias correction
% FORMAT P_trans = run_vbq_b1map(jobsubj)
%    jobsubj - are parameters for one subject out of the job list. 
%    NB: ONE SINGLE DATA SET FROM ONE SINGLE SUBJECT IS PROCESSED HERE, 
%    LOOP OVER SUBJECTS DONE AT HIGHER LEVEL. 
%    P_trans - a vector of file names with P_trans(1,:) = anatomical volume
%        for coregistration and P_trans(2,:) = B1 map in percent units.
%_______________________________________________________________________
% Written by E. Balteau, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
%_______________________________________________________________________
% $Id$

% retrieve b1_type from job and pass it as default value for the current
% processing:
b1_prot = jobsubj.b1_type;
vbq_get_defaults('b1map.b1_type.val',b1_prot);
% load the resulting default parameters:
b1map_defs = vbq_get_defaults(['b1map.',b1_prot,'.b1proc']);

% % force recalculation of default values taking b1_type into account:
% vbq_defaults;

P_trans = [];

if ~b1map_defs.avail
    disp('----- No B1 correction performed (semi-quantitative maps only) -----');
    return;
end

% calculate the B1 map if required
if b1map_defs.procreq
    if strfind(b1map_defs.b1_type.val,'AFI')
        % processing B1 map from AFI data
        P_trans  = calc_AFI_b1map(jobsubj, b1map_defs);
        
    elseif strfind(b1map_defs.b1_type.val,'EPI')
        % processing B1 map from SE/STE EPI data
        P_trans  = calc_SESTE_b1map(jobsubj, b1map_defs);
        
    end
else
    % return pre-processed B1 map
    disp('----- Loading pre-processed B1 map -----');
    P    = char(jobsubj.raw_fld.b1);   % the B1 map
    P_trans  = P(1:2,:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P_trans = calc_AFI_b1map(jobsubj, b1map_defs)

disp('----- Calculation of B1 map (AFI protocol) -----');

% NB: both phase and magnitude images can be provided but only the
% magnitude images (first series) are used. Phase images (second series)
% are not used. In each series, first image = TR2 (long TR) and second
% image = TR1 (short TR).
P = char(jobsubj.raw_fld.b1);   % 2 or 4 images
fileTR1 = P(2,:);
fileTR2 = P(1,:);
V1 = spm_vol(fileTR1);
V2 = spm_vol(fileTR2);
Vol1 = spm_read_vols(V1);
Vol2 = spm_read_vols(V2);

TR1 = 1; % only the ratio [TR2/TR1=n] matters
TR2 = b1map_defs.TR2TR1ratio;
alphanom = b1map_defs.alphanom;

% Mask = squeeze(Vol1);
% threshold = (prctile(Mask(:),98)-prctile(Mask(:),2))*0.1+prctile(Mask(:),2);
% Mask = (Mask>threshold);

B1map = acos((Vol2./Vol1*TR2/TR1-1)./(TR2/TR1*ones(size(Vol1))-Vol2./Vol1))*180/pi;
B1map_norm = abs(B1map)*100/alphanom;

% smoothed map
smB1map_norm = zeros(size(B1map_norm));
pxs = sqrt(sum(V1.mat(1:3,1:3).^2)); % Voxel resolution
smth = 8./pxs;
spm_smooth(B1map_norm,smB1map_norm,smth);

% masking
% B1map = B1map.*Mask;
% B1map_norm = B1map_norm.*Mask;
% smB1map_norm = smB1map_norm.*Mask;

%-Save everything in OUTPUT dir
%-----------------------------------------------------------------------
% determine output directory path
if isfield(jobsubj.output,'indir')
    outpath = fileparts(V1.fname);
else
    outpath = jobsubj.output.outdir{1};
end

[~, sname] = fileparts(V1.fname);

VB1 = V1;
% VB1.pinfo = [max(B1map(:))/16384;0;0];
% VB1.fname = fullfile(outpath, [sname '_B1map.nii']);
% spm_write_vol(VB1,B1map);

% VB1.pinfo = [max(B1map_norm(:))/16384;0;0];
% VB1.fname = fullfile(outpath, [sname '_B1map_norm.nii']);
% spm_write_vol(VB1,B1map_norm);

VB1.pinfo = [max(smB1map_norm(:))/16384;0;0];
VB1.fname = fullfile(outpath, [sname '_smB1map_norm.nii']);
spm_write_vol(VB1,smB1map_norm);

% requires anatomic image + map
P_trans  = char(char(fileTR1),char(VB1.fname));

% VB1.fname = fullfile(outpath, [sname '_B1map_mask.nii']);
% spm_write_vol(VB1,Mask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P_trans = calc_SESTE_b1map(jobsubj, b1map_defs)

disp('----- Calculation of B1 map (SE/STE EPI protocol) -----');

P    = char(jobsubj.raw_fld.b1); % B1 data - 11 pairs
Q    = char(jobsubj.raw_fld.b0); % B0 data - 3 volumes

V = spm_vol(P);
n = numel(V);
Y_tmptmp = zeros([V(1).dim(1:2) n]);
Y_ab = zeros(V(1).dim(1:3));
Y_cd = zeros(V(1).dim(1:3));
Index_Matrix = zeros([V(1).dim(1:3) b1map_defs.Nonominalvalues]);
real_Y_tmp = zeros([V(1).dim(1:2) 2*b1map_defs.Nonominalvalues]);
Y_tmp = zeros([V(1).dim(1:2) 2*b1map_defs.Nonominalvalues]);

Temp_matrix = zeros([V(1).dim(1:3) n/2]);
Ssq_matrix = zeros(V(1).dim(1:3));

for i=1:n/2
    Temp_matrix(:,:,:,i)=spm_read_vols(V(2*(i-1)+1));
end
Ssq_matrix=sqrt(sum(Temp_matrix.^2,4));

%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',V(1).dim(3),'B1 map fit','planes completed');

%-Loop over planes computing result Y
%-----------------------------------------------------------------------
clear Temp_mat;
corr_fact = exp(b1map_defs.TM/b1map_defs.T1);
for p = 1:V(1).dim(3),%loop over the partition dimension of the data set
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
    for i = 1:n/2
        M = inv(B*inv(V(1).mat)*V(1).mat);
        Y_tmptmp(:,:,((i-1)*2+1))  = real(acos(corr_fact*spm_slice_vol(V((i-1)*2+2),M,V(1).dim(1:2),0)./(spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0)+b1map_defs.eps))/pi*180/b1map_defs.beta(i)); % nearest neighbor interpolation
        Y_tmptmp(:,:,((i-1)*2+2))  = 180/b1map_defs.beta(i) - Y_tmptmp(:,:,((i-1)*2+1));
        Temp_mat(:,:,i) = spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0);
    end
    
    [~,indexes] = sort(Temp_mat,3);
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
            for k=1:b1map_defs.Nonominalvalues
                real_Y_tmp(x_nr,y_nr,2*k-1)=Y_tmptmp(x_nr,y_nr,2*indexes(x_nr,y_nr,n/2-k+1)-1);
                real_Y_tmp(x_nr,y_nr,2*k)=Y_tmptmp(x_nr,y_nr,2*indexes(x_nr,y_nr,n/2-k+1));
                Index_Matrix(x_nr,y_nr,p,k)=indexes(x_nr,y_nr,indexes(x_nr,y_nr,n/2-k+1));
            end
        end
    end
    
    Y_tmp = sort(real(real_Y_tmp), 3); % take the real value due to noise problems
    Y_sd  = zeros([V(1).dim(1:2) (b1map_defs.Nonominalvalues+1)]);
    Y_mn  = zeros([V(1).dim(1:2) (b1map_defs.Nonominalvalues+1)]);
    for i = 1:(b1map_defs.Nonominalvalues+1)
        Y_sd(:,:,i) = std(Y_tmp(:,:,i:(i + b1map_defs.Nonominalvalues-1)), [], 3);
        Y_mn(:,:,i) = mean(Y_tmp(:,:,i:(i + b1map_defs.Nonominalvalues-1)), 3);
    end
    
    [~,min_index] = min(Y_sd,[],3); % !! min_index is a 2D array. Size given by resolution along read and phase directions
    for x_nr = 1:V(1).dim(1)
        for y_nr = 1:V(1).dim(2)
            Y_ab(x_nr,y_nr,p) = Y_mn(x_nr,y_nr, min_index(x_nr,y_nr)); % Y_ab is the relative flip angle value averaged over the n flip angles (determined by minizing the SD i.e. keeping the most uniform relative flip angle values)
            Y_cd(x_nr,y_nr,p) = Y_sd(x_nr,y_nr, min_index(x_nr,y_nr)); % Y_cd is the corresponding standard deviation between the relative flip angle values
        end
    end
    spm_progress_bar('Set',p);
end

%-Save everything in OUTPUT dir
%-----------------------------------------------------------------------
% determine output directory path
if isfield(jobsubj.output,'indir')
    outpath = fileparts(V(1).fname);
else
    outpath = jobsubj.output.outdir{1};
end

V_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','B1 map [%]');
[~,name,e] = fileparts(V_save.fname);
P_B1 = fullfile(outpath,['B1map_' name e]);
V_save.fname = P_B1;
V_save = spm_write_vol(V_save,Y_ab*100);

W_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SD [%]');
P_SDB1 = fullfile(outpath,['SDmap_' name e]);
W_save.fname = P_SDB1;
W_save = spm_write_vol(W_save,Y_cd*100);

X_save = struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SE SSQ matrix');
P_SsqMat = fullfile(outpath,['SumOfSq' e]);
X_save.fname = P_SsqMat;
X_save = spm_write_vol(X_save,Ssq_matrix);


%-B0 undistortion
%-----------------------------------------------------------------------
% load default parameters and customize...
b0proc_defs = vbq_get_defaults('b0proc');
pm_defaults;
pm_defs = pm_def;
pm_defs.SHORT_ECHO_TIME = b0proc_defs.shorTE;
pm_defs.LONG_ECHO_TIME = b0proc_defs.longTE;
pm_defs.blipdir = b1map_defs.blipDIR;
pm_defs.TOTAL_EPI_READOUT_TIME = b1map_defs.EchoSpacing*b1map_defs.nPEacq;
pm_defs.MASKBRAIN = b0proc_defs.b0maskbrain;
pm_defs.STDTHRESH = 5;
pm_defs.HZTHRESH = b0proc_defs.HZTHRESH;
pm_defs.ERODEB1 = b0proc_defs.ERODEB1;
pm_defs.PADB1 = b0proc_defs.PADB1 ;
pm_defs.B1FWHM = b0proc_defs.B1FWHM; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
pm_defs.match_vdm = b0proc_defs.match_vdm;
pm_defs.pedir = b1map_defs.PEDIR;

mag1 = spm_vol(Q(1,:));
mag = mag1.fname;
phase1 = spm_vol(Q(3,:));
phase = phase1.fname;
scphase = FieldMap('Scale',phase);
fm_imgs = char(scphase.fname,mag);

anat_img1 = spm_vol(P_SsqMat);

[path,name,e] = fileparts(anat_img1.fname);
anat_img = {strcat(path,filesep,name,e)};
other_img{1} = char(V_save.fname);
other_img{2} = char(W_save.fname);

[fmap_img,unwarp_img] = B1Map_unwarp(fm_imgs,anat_img,other_img,pm_defs);
uanat_img{1}=unwarp_img{1}.fname;
ub1_img{1} = unwarp_img{2}.fname;
ustd_img{1} = unwarp_img{3}.fname;
fpm_img{1} = fmap_img{1};
vdm_img{1} = fmap_img{2};
[allub1_img] = B1Map_process(uanat_img,ub1_img,ustd_img,vdm_img,fpm_img,pm_defs);

P_trans  = char(char(uanat_img),char(allub1_img{2}.fname));

end