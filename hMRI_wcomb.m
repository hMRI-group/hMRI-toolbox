function hMRI_wcomb(PIn1,PIn2,Pw1,Pw2,PVG,PMSK)
% This function combines two input images (PIn1 and PIn2) using two weight 
% images (Pw1 and Pw2) for each input image, respectively. PVG is used as
% reference (space defining image). If PVG is unspecified, it will
% automatically use PIn1.
% S. Mohammadi 18/10/2019
% In:
% PIn           - Filepath & name of two input images
% Pw            - Filepath & name of two weight images
% 
% Out:
% 

% defaults -> here we need to the mpm_default settings
if(~exist('res','var'))
    res = 1;
end
if(~exist('kt','var'))
    kt = 20;
end
if(~exist('pval','var'))
    pval = 0.70;
end
if(~exist('dummy_am','var'))
    dummy_am = true;
end
if(~exist('smthk','var'))
    smthk = 1;
end
if(~exist('dim','var'))
    dim = 2;
end
dt = [spm_type('float32'),spm_platform('bigend')]; % for nifti output


% read in data
switch dim
    case 1
       dplane = [2 3];
    case 2
       dplane = [1 3];
    case 3
       dplane = [1 2];
end
VIn1 = spm_vol(PIn1);
VIn2 = spm_vol(PIn2);
Vw1 = spm_vol(Pw1);
Vw2 = spm_vol(Pw2);
if exist('PMSK','var') && ~isempty(PMSK)
    VMSK = spm_vol(PMSK);
    AMSK = hMRI_read_vols(VMSK,VIn1(1),res,[],dim);
else
    AMSK = ones(VIn1(1).dim);
end
for inx = 1:size(VIn1,1)
    % define reference volume
    if exist('PVG','var') && ~isempty(PVG)
        VG = spm_vol(PVG);
    else
        VG = spm_vol(VIn1(inx));    
    end
    % define output volume
    Pout = spm_file(VIn1(inx).fname,'prefix','wa_');
    
    Ntmp = hMRI_create_nifti(Pout,VG,dt,deblank([VIn1(inx).descrip  ' - weighted combination']));
    if dummy_am==true
        Pout = spm_file(VIn1(inx).fname,'prefix','am_');
        Ntmpam = hMRI_create_nifti(Pout,VG,dt,deblank([VIn1(inx).descrip  ' - arithmetic combination']));
    end
    spm_progress_bar('Init',VG.dim(3),Ntmp.descrip,'planes completed');
    
    % smooth weights
    vxg         = sqrt(sum(VG.mat(1:3,1:3).^2));
    smthk       = smthk.*vxg;
    Aw1 = hMRI_read_vols(Vw1(inx),VG,res,[],dim);
    Aw2 = hMRI_read_vols(Vw2(inx),VG,res,[],dim);
    sAw1 = Aw1;
    sAw2 = Aw2;
    spm_smooth (sAw1,Aw1, smthk);
    spm_smooth (sAw2,Aw2, smthk);
    
    for p = 1:VG.dim(dim)
        AIn1 = hMRI_read_vols(VIn1(inx),VG,res,p,dim);
        AIn2 = hMRI_read_vols(VIn2(inx),VG,res,p,dim);
%         Aw1 = hMRI_read_vols(Vw1(inx),VG,res,p);
%         Aw2 = hMRI_read_vols(Vw2(inx),VG,res,p);
        Aw1 = sAw1(:,:,p);
        Aw2 = sAw2(:,:,p);
        
        nAw1 = zeros(VG.dim(dplane));
        nAw2 = zeros(VG.dim(dplane));
        nAw1(abs(Aw1)>0) = abs(Aw1(abs(Aw1)>0))./(abs(Aw1(abs(Aw1)>0))+abs(Aw2(abs(Aw1)>0))+eps);        
        nAw2(abs(Aw2)>0) = abs(Aw2(abs(Aw2)>0))./(abs(Aw1(abs(Aw2)>0))+abs(Aw2(abs(Aw2)>0))+eps);

        f1 = local_fermi(nAw1(:),pval,kt,AMSK(:,:,p)); % we take 2- to down-weigh regions that have high res
        f2 = local_fermi(nAw2(:),pval,kt,AMSK(:,:,p)); % we take 2- to down-weigh regions that have high res
        Awavg = (AIn1(:).*f1 + AIn2(:).*f2)./(f1+f2);
        
        read_nifti_perm(Ntmp,reshape(Awavg,VG.dim(dplane)),dim,p)
        if dummy_am==true
            Aam = (AIn1(:) + AIn2(:))./2;
            read_nifti_perm(Ntmpam,reshape(Aam,VG.dim(dplane)),dim,p)    
        end
        
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');
end

end

function f= local_fermi(x,pval,kt,AMSK)
% This is a fermi function that goes from 0 to 2 at the point x0. kt
% provides the steepness with which the transition is done.
% S.Mohammadi 2.10.2019
    if isempty(find(AMSK>0))
        f = 1;
    else
        [ybin,xbin] = hist(x(AMSK>0),100);

        tmp = find(cumsum(ybin)>(1-pval)*sum(ybin));
        x0 = xbin(tmp(1));
        f= 2./(exp((abs(x)-x0).*kt)+1);
    end
end

function read_nifti_perm(Nif,A,dim,p)
% This function reads the data A into the nifti file Nif, accounting for
% the permutation of dimension defined by dim.
% S.Mohammadi 18.10.2019
    switch dim
        case 1
            Aout(1,:,:) = A; 
            Nif.dat(p,:,:) = Aout;
        case 2
            Aout(:,1,:) = A; 
            Nif.dat(:,p,:) = Aout;
        case 3
            Nif.dat(:,:,p) = A;
    end        
end