function PoutRO = hmri_wcomb_2mpms(PIn1,PIn2,Pw1,Pw2,outdir,PVG,PMSK)
% This function combines two input images (PIn1 and PIn2) using two weight
% images (Pw1 and Pw2) for each input image, respectively. PVG is used as
% reference (space defining image). If PVG is unspecified, it will
% automatically use PIn1.
% S. Mohammadi 18/10/2019
% In:
% PIn1   - Filepath & name of first input image
% PIn2   - Filepath & name of second input image
% Pw1    - Filepath & name of first weight image
% Pw2    - Filepath & name of second weight image
% outdir - Output directory
% PVG    - Filepath & name of reference image (optional: Image to which all
%          others will be resliced)
% PMSK   - Filepath & name of brain mask image

% Out:
% PoutRO - Robust combination

wcombparams = hmri_get_defaults('wcombparams');
write_am    = wcombparams.am; % write arithmetic mean image
write_error = wcombparams.errormaps; % write combined error maps
% Examples can be found in Mohammadi et al., NeuroImage, 2022, Supplementary Materials: S1: Efficiency of robust combination and the Fermi function
kt          = wcombparams.kt/100;
% The following parameters are for experts only and control the degree of
% smoothing of the errormaps perpendicular to dimension dim
smthk       = wcombparams.smthk;
dim         = wcombparams.dim;
res         = wcombparams.res;

dt = [spm_type('float32'),spm_platform('bigend')]; % for nifti output

% read in data
switch dim % choose which dimension slices will be read in; this is important if smoothing is used
    case 1
        dplane = [2 3];
    case 2
        dplane = [1 3];
    case 3
        dplane = [1 2];
end
VIn1 = spm_vol(PIn1);
VIn2 = spm_vol(PIn2);
Vw1  = spm_vol(Pw1);
Vw2  = spm_vol(Pw2);
% define reference volume
if exist('PVG','var') && ~isempty(PVG)
    VG = spm_vol(PVG);
else
    VG = spm_vol(VIn1);
end

if exist('PMSK','var') && ~isempty(PMSK)
    VMSK = spm_vol(PMSK);

    AMSK = false(VG.dim);
    for p = 1:VG.dim(dim)
        AMSK(:,:,p) = hmri_read_vols(VMSK,VG,p,res)>0;
    end
else
    AMSK = true(VG.dim);
end

% define output volumes
PoutRO = spm_file(VIn1.fname,'path',outdir,'suffix','_RO');
Ntmp = hmri_create_nifti(PoutRO,VG,dt,deblank([VIn1.descrip  ' - robust combination']));

if write_error
    Pout = spm_file(Vw1.fname,'path',outdir,'suffix',['_RO_k' num2str(kt)]);
    Ntmperror = hmri_create_nifti(Pout,VG,dt,deblank([Vw1.descrip  ' - robust combination error maps']));
end
if write_am
    Pout = spm_file(VIn1.fname,'path',outdir,'suffix','_AM');
    Ntmpam = hmri_create_nifti(Pout,VG,dt,deblank([VIn1.descrip  ' - arithmetic mean']));
end

spm_progress_bar('Init',VG.dim(3),Ntmp.descrip,'planes completed');

if smthk>0
    % smooth weights
    vxg         = sqrt(sum(VG.mat(1:3,1:3).^2));
    smthk       = smthk.*vxg;

    for p = 1:VG.dim(dim)
        Aw1 = hmri_read_vols(Vw1,VG,p,res);
        Aw2 = hmri_read_vols(Vw2,VG,p,res);
    end
    sAw1 = Aw1;
    sAw2 = Aw2;
    spm_smooth(sAw1,Aw1, smthk);
    spm_smooth(sAw2,Aw2, smthk);
end
for p = 1:VG.dim(dim)
    AIn1 = hmri_read_vols(VIn1,VG,p,res);
    AIn2 = hmri_read_vols(VIn2,VG,p,res);
    if write_am
        Aam = (AIn1 + AIn2)./2;
        write_nifti_slice(Ntmpam,reshape(Aam(:),VG.dim(dplane)),dim,p)
    end
    if smthk>0
        Aw1 = sAw1(:,:,p);
        Aw2 = sAw2(:,:,p);
    else
        Aw1 = hmri_read_vols(Vw1,VG,p,res);
        Aw2 = hmri_read_vols(Vw2,VG,p,res);
    end

    if nnz(AMSK(:,:,p))<1e1
        f1 = 1;
    else
        nAw1 = Aw1(AMSK(:,:,p))./Aw2(AMSK(:,:,p));
        nAw1(nAw1<0)=1;
        f1 = local_fermi(nAw1,kt); % we take 1- to downweight regions that have high res
    end

    Awavg = zeros(VG.dim(dplane));
    if nnz(AMSK(:,:,p))>0
        Awavg(AMSK(:,:,p)) = (AIn1(AMSK(:,:,p)).*f1 + AIn2(AMSK(:,:,p)).*(max(f1)-f1))./max(f1);
    end
    write_nifti_slice(Ntmp,reshape(Awavg,VG.dim(dplane)),dim,p)

    if write_error
        Awerr = zeros(VG.dim(dplane));
        Awerr(AMSK(:,:,p)) = (Aw1(AMSK(:,:,p)).*f1 + Aw2(AMSK(:,:,p)).*(max(f1)-f1))./max(f1);
        write_nifti_slice(Ntmperror,reshape(Awerr,VG.dim(dplane)),dim,p)
    end

    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

end

function f = local_fermi(x,kt)
% This is a fermi function that goes from 1 to 0 at the point 1. kt provides the steepness
% with which the transition is done.
% S.Mohammadi 2.10.2019
    f = 1./(exp((x-1)/kt)+1);
end

function write_nifti_slice(Nif,A,dim,p)
% This function writes the data A into the nifti file Nif, accounting for
% the choice of slice dimension defined by dim.
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