function [fR1, fR2s, fMT, fA, PPDw, PT1w]  = vbq_MTProtQA(P_mtw, P_pdw, P_t1w, TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w, fa_mtw, fa_pdw, fa_t1w, P_trans, P_receiv, flatAmap)

% function [fR1, fR2s, fMT, fA, PPDw, PT1w]  = vbq_MTProtQA(P_mtw, P_pdw, P_t1w, TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w, fa_mtw, fa_pdw, fa_t1w, P_trans, P_receiv, flatAmap)
% Evaluation function for multi-contrast multi-echo FLASH protocol
% P_mtw, P_pdw, P_t1w: MTw, PDw, T1w images (can be multiple echoes = images)
% TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w: echo times and TR of images
% fa_mtw, fa_pdw, fa_t1w: excitation flip angles of images
% P_trans: transmit bias map (p.u.) of B1 field of RF body coil
% P_receiv: sensitivity map of phased-array coil relative to BC (p.u.)
% flatAmap: default is 1 = flatten amplitude maps (0 == no flattening). Note that a B1 field map must be available for this feature.
%
% Gunther Helms, MR-Research in Neurology and Psychiatry, University of Goettingen
% Nikolaus Weiskopf, Antoine Lutti, John Ashburner, Wellcome Trust Centre for Neuroimaging at UCL, London


%
% Antoine Lutti 15/01/09
% This version of MTProt corrects for imperfect RF spoiling when a B1 map is loaded (line 229 and below)
% based on Preibisch and Deichmann's paper MRM 61:125-135 (2009).
% The values for P2_a and P2_b were obtained using the code supplied by
% Deichmann with the experimental parameters used to get our PDw and T1w
% images.
%

% $Id$

P_receiv = []; %CPL to supress error message

if nargin < 15
    flatAmap=[];
end
if nargin < 14
    P_receiv = [];
end
if nargin < 13
    P_trans = [];
end

if isempty(flatAmap) % by default flatten Amplitude maps
    flatAmap=0;
end

TE_limit = 30; % TE time up to which echoes are averaged (in ms)

% Settings for R.Deichmann steady state correction using T2=64ms at 3T
% Correction parameters were calculated for 3 different parameter sets:
if (feq(TR_pdw, 23.7) && feq(TR_t1w, 18.7) && feq(fa_pdw, 6) && feq(fa_t1w, 20))
    % 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
    disp('Classic FIL protocol');
    P2_a= [78.9228195006542,-101.113338489192,47.8783287525126];
    P2_b=[-0.147476233142129,0.126487385091045,0.956824374979504];
elseif (feq(TR_pdw, 24.5) && feq(TR_t1w, 24.5) && feq(fa_pdw, 5) && feq(fa_t1w, 29))
    % 2) new FIL/Helms protocol
    % PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
    disp('New FIL/Helms protocol');
    P2_a= [93.455034845930480,-120.5752858196904,55.911077913369060];
    P2_b=[-0.167301931434861,0.113507432776106,0.961765216743606];
elseif (feq(TR_pdw, 24.0) && feq(TR_t1w, 19.0) && feq(fa_pdw, 6) && feq(fa_t1w, 20))
    % 3) Siemens product sequence protocol used in Lausanne (G Krueger)
    %PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
    disp('Siemens product Lausanne (GK) protocol');
    P2_a= [67.023102027100880,-86.834117103841540,43.815818592349870];
    P2_b=[-0.130876849571103,0.117721807209409,0.959180058389875];
elseif (feq(TR_pdw, 23.7) && feq(TR_t1w, 23.7) && feq(fa_pdw, 6) && feq(fa_t1w, 28))
    % 4) High-res (0.8mm) FIL protocol:
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
    disp('High-res FIL protocol');
    P2_a= [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
    P2_b=[-0.218804328507184,0.178745853134922,0.939514554747592];
else
    warning('Warning!!! Spoiling correction not defined for this protocol. Defaulting to classic FIL protocol.');
    % 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
    P2_a= [78.9228195006542,-101.113338489192,47.8783287525126];
    P2_b=[-0.147476233142129,0.126487385091045,0.956824374979504];
end

% check that echo times are identical
for nr=1:min([length(TE_mtw), length(TE_pdw), length(TE_t1w)])
    if (TE_mtw(nr) == TE_pdw(nr)) & (TE_pdw(nr) == TE_t1w(nr))
        ok=1;
    else
        error('Echo times do not match!')
    end
end
nr_c_echoes = min([length(TE_mtw), length(TE_pdw), length(TE_t1w)]);

V_mtw   = spm_vol(P_mtw);
V_pdw   = spm_vol(P_pdw);
V_t1w   = spm_vol(P_t1w);

if ~isempty(P_trans) % B1 map (p.u.)
    V_trans   = spm_vol(P_trans);
else
    V_trans   = [];
end
if ~isempty(P_receiv) % sensitivity map (p.u.)
    V_receiv   = spm_vol(P_receiv);
else
    V_receiv   = [];
end

V_templ = spm_vol(P_mtw);
V       = V_templ(1);

[pth,nam,ext] = fileparts(P_mtw(1,:));

% calculate T2* map from PD echoes
dm        = V.dim;
spm_progress_bar('Init',dm(3),'R2* fit','planes completed');
disp('----- Calculation of T2* map -----');
V         = V_templ(1);
n         = numel(V_pdw);
dt        = [spm_type('float32'),spm_platform('bigend')];
Ni        = nifti;
Ni.mat    = V.mat;
Ni.mat0   = V.mat;
Ni.descrip='R2* map [1/ms]';
Ni.dat    = file_array(fullfile(pth,[nam '_R2s' '.nii']),dm,dt, 0,1,0);
create(Ni);
fR2s = fullfile(pth,[nam '_R2s' '.nii']);

reg = [ones(numel(TE_pdw),1) TE_pdw(:)];
W   = (reg'*reg)\reg';
W   = -W(2,:)';
for p = 1:dm(3),
    M = spm_matrix([0 0 p 0 0 0 1 1 1]);
    Y = zeros(dm(1:2));
    for i = 1:n,
        M1 = V_pdw(i).mat\V_pdw(1).mat*M;
        Y  = Y + W(i)*log(max(spm_slice_vol(V_pdw(i),M1,dm(1:2),1),1));
    end
    Ni.dat(:,:,p) = max(min(Y,10),-10); % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

% Average first few echoes for increased SNR and fit T2*
disp('----- Reading and averaging the images -----');

nr_TE_limit = min(find(TE_mtw > TE_limit));
avg_nr      = min([nr_c_echoes nr_TE_limit]);

PP   = {P_mtw,P_pdw,P_t1w};
nam1 = {'MTw','PDw','T1w'};
avg  = [0 0 0];
for ii=1:3,
    V         = spm_vol(PP{ii});
    dm        = V(1).dim;
    Ni        = nifti;
    Ni.mat    = V(1).mat;
    Ni.mat0   = V(1).mat;
    Ni.descrip= sprintf('Averaged %s images', nam1{ii});
    Ni.dat    = file_array(fullfile(pth,[nam '_' nam1{ii} '.nii']),dm,dt, 0,1,0);
    create(Ni);
    spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
    sm = 0;
    for p=1:dm(3),
        M = spm_matrix([0 0 p]);
        Y = zeros(dm(1:2));
        for nr=1:avg_nr,
            M1 = V(nr).mat\V(1).mat*M;
            Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),1);
        end
        Ni.dat(:,:,p) = Y/avg_nr;
        sm = sm + sum(Y(:))/avg_nr;
        spm_progress_bar('Set',p);
    end
    avg(ii) = sm/prod(dm);
    spm_progress_bar('Clear');
end
fPD = fullfile(pth,[nam '_' nam1{ii} '.nii']);

PPDw = fullfile(pth,[nam '_PDw' '.nii']);
PMTw = fullfile(pth,[nam '_MTw' '.nii']);
PT1w = fullfile(pth,[nam '_T1w' '.nii']);

if true,
    disp('----- Coregistering the images -----');
    coreg_mt(PPDw, PMTw);
    coreg_mt(PPDw, PT1w);
    if ~isempty(V_trans)
        coreg_bias_map(PPDw, P_trans);
    end
    if ~isempty(V_receiv)
        coreg_bias_map(PPDw, P_receiv);
    end
end

VPDw = spm_vol(PPDw);
VMTw = spm_vol(PMTw);
VT1w = spm_vol(PT1w);
if ~isempty(V_trans)
    Vtrans = spm_vol(P_trans(2,:)); % map only, we do not need the T1w image any more
else
    Vtrans = [];
end
if ~isempty(V_receiv)
    Vreceiv = spm_vol(P_receiv(2,:)); % map only, we do not need the T1w image any more
else
    Vreceiv = [];
end


nam2    = {'R1','A','MT','MTR_synt'};
descrip = {'R1 map [1000/s]', 'A map','Delta MT map', 'Synthetic MTR image'};
if (TR_mtw == TR_pdw) & (fa_mtw == fa_pdw),
    nam2    = {nam2{:}, 'MTR','MTRdiff'};
    descrip = {descrip{:}, 'Classic MTR image','Percent diff. MTR image (RD/BD)'};
end
if (flatAmap)
    nam2    = {nam2{:}, 'MTforA'};
    descrip = {descrip{:}, 'B1-uncorrected MT map for A flattening'};
end
Nmap    = nifti;
for ii=1:numel(nam2),
    V         = V_templ(1);
    dm        = V(1).dim;
    Ni        = nifti;
    Ni.mat    = V(1).mat;
    Ni.mat0   = V(1).mat;
    Ni.descrip= descrip{ii};
    Ni.dat    = file_array(fullfile(pth,[nam '_' nam2{ii} '.nii']),dm,dt, 0,1,0);
    create(Ni);
    Nmap(ii) = Ni;
end
fR1 = fullfile(pth,[nam '_' nam2{1} '.nii']);
fA  = fullfile(pth,[nam '_' nam2{2} '.nii']);
fMT = fullfile(pth,[nam '_' nam2{3} '.nii']);

disp('----- Calculating the maps -----');
M0 = Ni.mat;
dm = size(Ni.dat);

fa_pdw = fa_pdw * pi / 180;
fa_mtw = fa_mtw * pi / 180;
fa_t1w = fa_t1w * pi / 180;

spm_progress_bar('Init',dm(3),'Calculating maps','planes completed');

for p=1:dm(3),
    M = M0*spm_matrix([0 0 p]);

    MTw = spm_slice_vol(VMTw,VMTw.mat\M,dm(1:2),3);
    PDw = spm_slice_vol(VPDw,VPDw.mat\M,dm(1:2),3);
    T1w = spm_slice_vol(VT1w,VT1w.mat\M,dm(1:2),3);
    
    if ~isempty(Vtrans)
        f_T = spm_slice_vol(Vtrans,Vtrans.mat\M,dm(1:2),3)/100; % divide by 100, since p.u. maps
    else
        f_T = [];
    end
    if ~isempty(Vreceiv)&~isempty(Vtrans)
        f_R = spm_slice_vol(Vreceiv,Vreceiv.mat\M,dm(1:2),3)/100; % divide by 100, since p.u. maps
        f_R = f_R .* f_T; % f_R is only the sensitivity map and not the true receive bias map, therefore needs to be multiplied by transmit bias (B1+ approx. B1- map)
    else
        f_R = [];
    end
    
    % Standard magnetization transfer ratio (MTR) in percent units [p.u.]
    % only if  trpd = trmt and fapd = fmt
    % else calculate "synthetic MTR using A and T1 (see below)
    if numel(Nmap)>4&&(TR_mtw == TR_pdw) && (fa_mtw == fa_pdw),
        MTR = (PDw-MTw)./(PDw+eps) * 100;
        % write MTR image
        Nmap(5).dat(:,:,p) = max(min(MTR,50),-50);
        
        % calculate a modified MTR map according to RD/BD
        MTR = 100*(PDw-MTw)./(eps+PDw).*(MTw./(eps+PDw)<1.3&MTw./(eps+PDw)>0&PDw>25);
        Nmap(6).dat(:,:,p) = max(min(MTR,50),-50);
    end
    
    % calculating T1 and A from a rational approximation of the Ernst equation using radian units
    % divide by fa and subtract
    % PD_d = PDw / fa_pdw;
    % T1_d = T1w / fa_t1w;
    %PD_T1_d = (PDw / fa_pdw) - (T1w / fa_t1w);
    
    % multiply by fa and divide by 2TR and subtract
    % PD_m = PDw * fa_pdw / 2 / TR_pdw;
    % T1_m = T1w * fa_t1w / 2 / TR_t1w;   % nw: corrected from T1_d to T1_m, correct?!,
    %T1_PD_m = (T1w * fa_t1w / 2 / TR_t1w) - (PDw * fa_pdw / 2 / TR_pdw);
    
    
    if isempty(f_T)
        % semi-quantitative T1
        T1 = ((PDw / fa_pdw) - (T1w / fa_t1w)) ./ ...
            max((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw),eps);
        R1 = (((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw)) ./ ...
            max(((PDw / fa_pdw) - (T1w / fa_t1w)),eps))*10^6;
    else
        % Transmit bias corrected quantitative T1 values
        % correct T1 for transmit bias f_T with fa_true = f_T * fa_nom
        % T1corr = T1 / f_T / f_T
        T1 = P2_a(1)*f_T.^2+P2_a(2)*f_T+P2_a(3)+(P2_b(1)*f_T.^2+P2_b(2)*f_T+P2_b(3)).*((((PDw / fa_pdw) - (T1w / fa_t1w)+eps) ./ ...
            max((T1w * fa_t1w / 2 / TR_t1w) - (PDw * fa_pdw / 2 / TR_pdw),eps))./f_T.^2);
        
        R1=1./T1*10^6;
        
        %          R1App=f_T.^2.*(((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw)) ./ ...
        %              max(((PDw / fa_pdw) - (T1w / fa_t1w)),eps));
        %          A_poly=P2_a(1)*f_T.^2+P2_a(2)*f_T+P2_a(3);B_poly=(P2_b(1)*f_T.^2+P2_b(2)*f_T+P2_b(3));
        %          R1 = R1App./(max(R1App.*A_poly+B_poly,eps))*10^6;
    end
    T1_forMT = ((PDw / fa_pdw) - (T1w / fa_t1w)) ./ ...
        max((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw),eps);
    T1       = max(T1,0);

    
    R1(R1<0)=0;
    tmp      = R1;
    Nmap(1).dat(:,:,p) = min(max(tmp,-10000),10000); %Truncating images to dynamic range. Originally 20000. Negative values are allowed or min(max(tmp,0),2000)
    
    % A values proportional to PD
    if (~isempty(f_T)) && (~isempty(f_R))
        % Transmit and receive bias corrected quantitative A values
        % again: correct A for transmit bias f_T and receive bias f_R
        % Acorr = A / f_T / f_R , proportional PD
        A = (T1 .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w))./f_T./f_R;
    elseif(~isempty(f_T))&&(isempty(f_R))&&flatAmap
        A = T1 .* (T1w .*(fa_t1w*f_T) / 2 / TR_t1w) + (T1w ./ (fa_t1w*f_T));
    else
        % semi-quantitative A
        A = T1 .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w);
    end
    
    A_forMT = T1_forMT .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w);
    tmp      = A;
    Athresh=10^5;
    Nmap(2).dat(:,:,p) = max(min(tmp,Athresh),-Athresh); % dynamic range increased to 10^5 to accommodate phased-array coils and symmetrical for noise distribution
    
    % MT in [p.u.]; offset by - famt * famt / 2 * 100 where MT_w = 0 (outside mask)
    MT       = ( (A_forMT * fa_mtw - MTw) ./ (MTw+eps) ./ (T1_forMT + eps) * TR_mtw - fa_mtw * fa_mtw / 2 ) * 100;
    if (~isempty(f_T))
        MT = MT .* (1 - 0.4) ./ (1 - 0.4 * f_T);
    end
    if flatAmap%The two outer voxels in all directions are nulled in order to remove artefactual effects from the MT map
        MTforA=MT;
        if (p==1)||(p==2)||(p==dm(3))||(p==dm(3)-1)
            MTforA=zeros(size(MT,1),size(MT,2));
        else
            MTforA(1,:)=0;MTforA(2,:)=0;MTforA(end-1,:)=0;MTforA(end,:)=0;
            MTforA(:,1)=0;MTforA(:,2)=0;MTforA(:,end-1)=0;MTforA(:,end)=0;
        end
        tmp      = MTforA;
        if (TR_mtw == TR_pdw) & (fa_mtw == fa_pdw)
            Nmap(7).dat(:,:,p) = max(min(tmp,50),-50);
        else
            Nmap(5).dat(:,:,p) = max(min(tmp,50),-50);
        end
    end
    
    tmp      = MT;
    Nmap(3).dat(:,:,p) = max(min(tmp,5),-5);
    
    % calculate synthetic reference signal at trmt and famt using the
    % rational approximation of the Ernst equation
    S_ref      = A_forMT .* fa_mtw * TR_mtw ./ (T1_forMT+eps) ./ ( TR_mtw ./ (T1_forMT+eps) +  fa_mtw * fa_mtw / 2 );
    % MTR_synt = (S_ref ./ MTw - 1) * 100;
    MTR_synt   = (S_ref-MTw) ./ (S_ref+eps) * 100;
    tmp      = MTR_synt;
    Nmap(4).dat(:,:,p) = max(min(tmp,50),-50);
    spm_progress_bar('Set',p);
    
        
end
if(~isempty(f_T))&&(isempty(f_R))&&flatAmap
    Aflattening(pth,Athresh)
end
spm_progress_bar('Clear');

return;

% --------------------------------------------------------------
function [] = coreg_mt(P_ref, P_src)
% coregisters the structural images
% for MT protocol

for src_nr=1:size(P_src, 1)
    P_src(src_nr,:);
    VG = spm_vol(P_ref);
    VF = spm_vol(P_src(src_nr,:));
    %coregflags.sep = [2 1];
    coregflags.sep = [4 2];
%     x = spm_coreg(VG,VF, coregflags);
    x = zeros(1,6);
    %x  = spm_coreg(mireg(i).VG, mireg(i).VF,flags.estimate);
    M  = inv(spm_matrix(x));
    MM = spm_get_space(deblank(VF.fname));
    spm_get_space(deblank(deblank(VF.fname)), M*MM);
end
return;


% -----------------------------------------------------------------

% --------------------------------------------------------------
function [] = coreg_bias_map(P_ref, P_src)
% coregisters the B1 or receive maps with
% the structurals in the MT protocol

P_src(1,:);
VG = spm_vol(P_ref);
VF = spm_vol(P_src(1,:));
%coregflags.sep = [2 1];
coregflags.sep = [4 2];
%x = spm_coreg(VG,VF, coregflags);
x = zeros(1,6);
%x  = spm_coreg(mireg(i).VG, mireg(i).VF,flags.estimate);
M  = inv(spm_matrix(x));
MM = spm_get_space(deblank(VF.fname));
spm_get_space(deblank(deblank(VF.fname)), M*MM);

VF2 = spm_vol(P_src(2,:)); % now also apply transform to the map
M  = inv(spm_matrix(x));
MM = spm_get_space(deblank(VF2.fname));
spm_get_space(deblank(deblank(VF2.fname)), M*MM);

return;

function Aflattening(pth,Athresh)
% TempT1=spm_select('FPList',pth ,'^.*_R1.(img|nii)$');
% myT1=10^6./spm_read_vols(spm_vol(TempT1));
% myT1(myT1>=10000)=0;myT1(isnan(myT1))=0;
% Corr_save=spm_vol(TempT1);
% [pth,name,e,v] = fileparts(TempT1);
% Corr_save.fname = fullfile(pth,['TempT1' e]);
% spm_write_vol(Corr_save,myT1);%Used further down for flattening of the A maps.


pm_defaults
myflags.template = pm_def.MFLAGS.TEMPLATE;
myflags.fwhm=pm_def.MFLAGS.FWHM;
myflags.nerode=pm_def.MFLAGS.NERODE;
myflags.ndilate=pm_def.MFLAGS.NDILATE;
myflags.thresh=pm_def.MFLAGS.THRESH;
myflags.reg = pm_def.MFLAGS.REG;
myflags.graphics=pm_def.MFLAGS.GRAPHICS;

myflags.fwhm=8;
myflags.nerode=2;
myflags.ndilate=6;
myflags.thresh=0.5;

T1formask=spm_select('FPList',pth ,'^.*_MTforA.(img|nii)$');
Corr_save=spm_vol(T1formask);
[pth,name,e] = fileparts(T1formask);
[Mask,Calibmask]=MaskT1(spm_vol(T1formask),myflags);
temp=spm_select('FPList',pth,'^(c|Bias).*\_MT');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end
delete(spm_select('FPList',pth,'^.*_MTforA_seg8'));
delete(spm_select('FPList',pth,'^.*_MTforA.(img|nii)$'));

P=spm_select('FPList',pth ,'^.*_A.(img|nii)$');
A=spm_read_vols(spm_vol(P));
A(A==Inf)=0;A(isnan(A))=0;A(A==Athresh)=0;

maskedA=A.*Mask;
Corr_save.fname = fullfile(pth,['maskedAmap' e]);
spm_write_vol(Corr_save,maskedA);%

maskedA_cell=cellstr(spm_select('FPList',pth,'^maskedAmap.*\.(img|nii)$'));
clear matlabbatch
matlabbatch{1}.spm.tools.preproc8.channel.vols = {maskedA_cell{1}};
matlabbatch{1}.spm.tools.preproc8.channel.write = [1 0];
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 10^(-5);
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 50;
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

% read bias field into memory
% -------------------------------------------------------------------------
bf = fullfile(pth, spm_select('List', pth, '^BiasField.*\.(img|nii)$'));
bfv = spm_vol(bf);
BF = double(spm_read_vols(bfv));

% apply bias field
% -------------------------------------------------------------------------
% read file
% fn = [P{1}];
V = spm_vol(P);
Y = spm_read_vols(V);
% apply bias field
Y = BF.*Y;
Calib_A=Calibmask.*Y;
% Y=Y/mean(Calib_A(Calib_A~=0))*100;%Calibmask is a CSF mask. Calibration of PD values based on the mean A value in CSF for each individual subject
% Y=Y/(7.25e3)*100;%Calibration of proton density values based on average
% CSF A values calculated on a bunch of subjects
Y=Y/mean(Calib_A(Calib_A~=0))*69;%Calibmask is a WM mask. Calibration on white matter from the litterature
sprintf('mean CSF: %04d',mean(Calib_A(Calib_A~=0)))
sprintf('SDs CSF %04d',std(Calib_A(Calib_A~=0),[],1))
Y(Y>200)=0;

% Y(Y>15000)=15000;
V.fname = P;
spm_write_vol(V,Y);

temp=spm_select('FPList',pth,'^.*\maskedAmap');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end


return;

function [bmask,CSF] = MaskT1(P,flags)

if nargin < 2 || isempty(flags)
    flags.template=fullfile(spm('Dir'),'templates','T1.nii');
    flags.fwhm=5;
    flags.nerode=2;
    flags.ndilate=4;
    flags.thresh=0.5;
    flags.reg = 0.02;
    flags.graphics=0;
end

disp('Segmenting and extracting brain...');
seg_flags.estimate.reg=flags.reg;
seg_flags.graphics = flags.graphics;


% Antoine 01/2011; Use new segment
epi_all=cellstr(P.fname);
% epi_all=cellstr(spm_select('FPList',pth,'^maskedAmap.*\.(img|nii)$'));
clear matlabbatch
matlabbatch{1}.spm.tools.preproc8.channel.vols = {epi_all{1}};
matlabbatch{1}.spm.tools.preproc8.channel.write = [1 0];
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
[pth,name,e] = fileparts(P.fname);
VO=spm_vol(cellstr(spm_select('FPList',pth,'^c.*\.(img|nii)$')));
GM=spm_read_vols(VO{1,1});WM=spm_read_vols(VO{2,1});CSF=spm_read_vols(VO{3,1});Other1=spm_read_vols(VO{4,1});Other2=spm_read_vols(VO{5,1});
bmask=zeros(size(GM,1),size(GM,2),size(GM,3));
bmask((GM+WM+CSF+Other1+Other2)>=flags.thresh)=1;
% CSF(CSF<0.95)=0;CSF(CSF>=0.95)=1;
CSF(WM<0.95)=0;CSF(WM>=0.95)=1;%CSF is actually a white matter mask
return;


%__________________________________________________________________________


function bl = feq(val, comp_val)
% floating point comparison
bl = abs(val - comp_val) <= eps(comp_val);
return;
