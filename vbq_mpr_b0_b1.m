function out=vbq_mpr_b0_b1(job)
% Calculation of B1 mapping data 3D EPI spin echo (SE) and stimulated (STE) echo images (see Jiru and Klose MRM 2006).
% Corresponding scanning protocol/sequence: al_B1mapping_v2a and
% al_B1mapping_v2b
% Input: 11 pairs of (SE, STE) images for B1 map calculation and 3 images
% for B0 map calculation.
% This macro calls the functions B1Map_unwarp and B1Map_process for
% correction of image distortions, padding and smoothing of the images.
% Output:
%     - distorted B1 (B1map_----) and error (SDmap_----) maps
%     - undistorted B1 (uB1map_----) and error (uSDmap_----) maps
%     - undistorted, masked and padded B1 maps (muB1map_---------)
%     - undistorted, masked, padded and smoothed B1 maps (smuB1map_---------) i.e. FULLY PROCESSED
% At each voxel, this macro selects the 5 pairs of (SE,STE image) (out of
% 11) with maximum signal amplitude in the SE images.
% The sum of square image of all SE images is created (SumOfSq) and
% undistorted (uSumOfSq) for coregistration of the B1 map to an anatomical dataset
% former B1map_v2.m

% $Id$

job=vbq_process_data_spec(job);

out.R1 = {};
out.R2s = {};
out.A = {};
out.MT = {};
out.T1w = {};

for in=1:numel(job.subj)
    local_job.subj=job.subj(in);
    out_temp=vbq_mpr_b0_b1_local(local_job,job.b1_type);
    out.subj(in)=out_temp.subj(1);
    out.R1{end+1} = out.subj(in).R1{1};
    out.R2s{end+1} = out.subj(in).R2s{1};
    out.MT{end+1} = out.subj(in).MT{1};
    out.A{end+1} = out.subj(in).A{1};
    out.T1w{end+1} = out.subj(in).T1w{1};
end
end

function out_loc = vbq_mpr_b0_b1_local(job,b1_type)

% b1_type = '3D_EPI_v2b'; 
% B1map_type = '3D_EPI_rect700';

if strcmp(b1_type, '3D_EPI_v2b')
    beta=115:-5:65;
    TM = 31.2;
    Nonominalvalues=5;
elseif strcmp(b1_type, '3D_EPI_v2b_long')
    beta=135:-5:65;
    TM = 33.24;
    Nonominalvalues=5;
elseif strcmp(b1_type, '3D_EPI_rect700')
    beta=80:5:100;
    TM = 33.53;
    Nonominalvalues=5;
elseif strcmp(b1_type, '3D_EPI_12ch')
    beta = 115:-7.5:62.5;
    TM = 31.2;
    Nonominalvalues = 3;
else
    warning('Warning!!! Unknown type of B1 map data. Defaulting to 3D_EPI_v2b');
    beta=115:-5:65;
    TM = 31.2;
    Nonominalvalues=5;
end
T1 = 1192; %ms, strictly valid only at 3T
eps=0.0001;


for ip=1:numel(job.subj)
    P    = char(job.subj(ip).raw_fld.b1);   % the B1 - 11 pairs
    Q    = char(job.subj(ip).raw_fld.b0);    % the B0 - 3 pairs
    
    % if nargin < 2
    %     P = spm_select(Inf,'image','Select images for B1 maps');
    %     Q = spm_select(Inf,'image','Select images for B0 map');
    % end
    disp('----- Calculation of B1 map -----');
    
    V=spm_vol(P);
    n   = numel(V);
    Y_tmptmp   = zeros([V(1).dim(1:2) n]);
    Y_ab = zeros(V(1).dim(1:3));
    Y_cd = zeros(V(1).dim(1:3));
    indexes=zeros(n/2,1);
    Index_Matrix = zeros([V(1).dim(1:3) Nonominalvalues]);
    real_Y_tmp   = zeros([V(1).dim(1:2) 2*Nonominalvalues]);
    Y_tmp   = zeros([V(1).dim(1:2) 2*Nonominalvalues]);
    
    Temp_matrix   = zeros([V(1).dim(1:3) n/2]);
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
    corr_fact = exp(TM/T1);
    for p = 1:V(1).dim(3),%loop over the partition dimension of the data set
        B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
        
        for i = 1:n/2
            M = inv(B*inv(V(1).mat)*V(1).mat);
            Y_tmptmp(:,:,((i-1)*2+1))  = real(acos(corr_fact*spm_slice_vol(V((i-1)*2+2),M,V(1).dim(1:2),0)./(spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0)+eps))/pi*180/beta(i)); % nearest neighbor interpolation
            Y_tmptmp(:,:,((i-1)*2+2))  = 180/beta(i) - Y_tmptmp(:,:,((i-1)*2+1));
            Temp_mat(:,:,i)=spm_slice_vol(V((i-1)*2+1),M,V(1).dim(1:2),0);
        end
        [~,indexes]=sort(Temp_mat,3);
        for x_nr = 1:V(1).dim(1)
            for y_nr = 1:V(1).dim(2)
                for k=1:Nonominalvalues
                    real_Y_tmp(x_nr,y_nr,2*k-1)=Y_tmptmp(x_nr,y_nr,2*indexes(x_nr,y_nr,n/2-k+1)-1);
                    real_Y_tmp(x_nr,y_nr,2*k)=Y_tmptmp(x_nr,y_nr,2*indexes(x_nr,y_nr,n/2-k+1));
                    Index_Matrix(x_nr,y_nr,p,k)=indexes(x_nr,y_nr,indexes(x_nr,y_nr,n/2-k+1));
                end
            end
        end
        Y_tmp = sort(real(real_Y_tmp), 3); % take the real value due to noise problems
        Y_sd   = zeros([V(1).dim(1:2) (Nonominalvalues+1)]);
        Y_mn   = zeros([V(1).dim(1:2) (Nonominalvalues+1)]);
        for i = 1:(Nonominalvalues+1)
            Y_sd(:,:,i) = std(Y_tmp(:,:,i:(i+Nonominalvalues-1)), [], 3);
            Y_mn(:,:,i) = mean(Y_tmp(:,:,i:(i+Nonominalvalues-1)), 3);
        end
        [~,min_index] = min(Y_sd,[],3);%!!min_index is a 2D array. Size given by resolution along read and phase directions
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
    
    % save workspace
    X_save=struct('fname',V(1).fname,'dim',V(1).dim,'mat',V(1).mat,'dt',V(1).dt,'descrip','SE SSQ matrix');
    [p,name,e] = fileparts(X_save.fname);
    P_SsqMat = fullfile(p,['SumOfSq' e]);
    X_save.fname = P_SsqMat;
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
    
    anat_img1=spm_vol(P_SsqMat);
    
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
    
    
    P_mtw    = char(job.subj(ip).raw_mpm.MT);
    P_pdw    = char(job.subj(ip).raw_mpm.PD);
    P_t1w    = char(job.subj(ip).raw_mpm.T1);
    P_trans  = strvcat(char(uanat_img),char(allub1_img{2}.fname));
    
    P_receiv = [];
    
    
    p = hinfo(P_mtw);
    TE_mtw = cat(1,p.te);
    TR_mtw = p(1).tr;
    fa_mtw = p(1).fa;
    
    p = hinfo(P_pdw);
    TE_pdw = cat(1,p.te);
    TR_pdw = p(1).tr;
    fa_pdw = p(1).fa;
    
    p = hinfo(P_t1w);
    TE_t1w = cat(1,p.te);
    TR_t1w = p(1).tr;
    fa_t1w = p(1).fa;
    
    [fR1, fR2s, fMT, fA, PPDw, PT1w]  = MTProt(P_mtw, P_pdw, P_t1w, TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w, fa_mtw, fa_pdw, fa_t1w, P_trans, P_receiv);
    if isfield(job.subj(ip).output,'indir') && job.subj(ip).output.indir == 1
        cwd = fileparts(fR1); %MFC - was set as p which is defined as hinfo(P_t1w) => use last path instead
    else
        % MFC: Only move files if a different directory is chosen - can't
        % move a file to itself...
        cwd=job.subj(ip).output.outdir{1};
        movefile(fR1,cwd);
        movefile(fR2s,cwd);
        movefile(fMT,cwd);
        movefile(fA,cwd);
        movefile(PPDw,cwd);
        movefile(PT1w,cwd);
    end
    
    out_loc.subj(ip).R1={fullfile(cwd,spm_str_manip(fR1,'t'))};
    out_loc.subj(ip).R2s={fullfile(cwd,spm_str_manip(fR2s,'t'))};
    out_loc.subj(ip).MT={fullfile(cwd,spm_str_manip(fMT,'t'))};
    out_loc.subj(ip).A={fullfile(cwd,spm_str_manip(fA,'t'))};
    out_loc.subj(ip).T1w={fullfile(cwd,spm_str_manip(PT1w,'t'))};
    
    f = fopen(fullfile(cwd, '_finished_'), 'wb');
    fclose(f);
end
end

function p = hinfo(P)
N = nifti(P);
for ii = 1:numel(N),
    tmp = regexp(N(ii).descrip,...
        'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg',...
        'names');
    p(ii).tr=str2num(tmp.tr);
    p(ii).te=str2num(tmp.te);
    p(ii).fa=str2num(tmp.fa);
end

end
