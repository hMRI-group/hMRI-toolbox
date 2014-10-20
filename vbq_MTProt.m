function [fR1, fR2s, fMT, fA, PPDw, PT1w]  = vbq_MTProt(P_mtw, P_pdw, P_t1w, P_trans, P_receiv)

% Evaluation function for multi-contrast multi-echo FLASH protocol
% P_mtw, P_pdw, P_t1w: MTw, PDw, T1w images (can be multiple echoes = images)
% TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w: echo times and TR of images
% fa_mtw, fa_pdw, fa_t1w: excitation flip angles of images
% P_trans: transmit bias map (p.u.) of B1 field of RF body coil
% P_receiv: sensitivity map of phased-array coil relative to BC (p.u.)
%
% Gunther Helms, MR-Research in Neurology and Psychiatry, University of Goettingen
% Nikolaus Weiskopf, Antoine Lutti, John Ashburner, Wellcome Trust Centre for Neuroimaging at UCL, London



%
% Antoine Lutti 15/01/09
% This version of MTProt corrects for imperfect RF spoiling when a B1 map
% is loaded (line 229 and below) 
% based on Preibisch and Deichmann's paper MRM 61:125-135 (2009).
% The values for P2_a and P2_b were obtained using the code supplied by
% Deichmann with the experimental parameters used to get our PDw and T1w
% images.
%

% MFC 31.05.2013    If the spoiling correction has not been defined for the
%                   protocol then none is applied.
%
% MFC 27.03.2014    Use PDw image as template for writing out maps to be
%                   robust in cases of inconsistent FoV positioning.

% $Id$

disp('----- Create maps from multi-contrast multi-echo FLASH protocol -----');

P_receiv = []; %CPL to supress error message

if nargin < 5
    P_receiv = [];
end
if nargin < 4
    P_trans = [];
end

TE_limit = 30; % TE time up to which echoes are averaged (in ms)


%% retrieves acquisition parameters and update defaults parameters set
p = hinfo(P_mtw);
vbq_get_defaults('MPMacq.TE_mtw', cat(1,p.te));
vbq_get_defaults('MPMacq.TR_mtw', p(1).tr);
vbq_get_defaults('MPMacq.fa_mtw', p(1).fa);

p = hinfo(P_pdw);
vbq_get_defaults('MPMacq.TE_pdw', cat(1,p.te));
vbq_get_defaults('MPMacq.TR_pdw', p(1).tr);
vbq_get_defaults('MPMacq.fa_pdw', p(1).fa);

p = hinfo(P_t1w);
vbq_get_defaults('MPMacq.TE_t1w', cat(1,p.te));
vbq_get_defaults('MPMacq.TR_t1w', p(1).tr);
vbq_get_defaults('MPMacq.fa_t1w', p(1).fa);

% force running default script to update steady-state correction parameters
vbq_defaults;


%% locally retrieves default parameters from vbq_defaults
% load threshold to save qMRI maps
threshall = vbq_get_defaults('qMRI_maps_thresh');
% load PD maps processing parameters
PDproc = vbq_get_defaults('PDproc');
% retrieve acquisition parameters
MPMacq = vbq_get_defaults('MPMacq');
% NB: for better readability, avoiding lengthy notations, the following
% parameters are re-written out of the MPMacq structure :\...
TE_pdw = MPMacq.TE_pdw; 
TE_mtw = MPMacq.TE_mtw; 
TE_t1w = MPMacq.TE_t1w; 
TR_pdw = MPMacq.TR_pdw; 
TR_mtw = MPMacq.TR_mtw; 
TR_t1w = MPMacq.TR_t1w; 
fa_pdw = MPMacq.fa_pdw; 
fa_mtw = MPMacq.fa_mtw; 
fa_t1w = MPMacq.fa_t1w; 
% RF spoiling correction parameters
RFC = vbq_get_defaults('rfcorr');

% a few words to summarize the situation...
disp(['INFO: Acquisition protocol = ' RFC.tag]);

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

V_templ = spm_vol(P_pdw);
V       = V_templ(1);

[pth,nam,ext] = fileparts(P_mtw(1,:));


%% calculate T2* map from PD echoes
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
    Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

% Average first few echoes for increased SNR and fit T2*
disp('----- Reading and averaging the images -----');

nr_TE_limit = find(TE_mtw > TE_limit,1);
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
if (PDproc.PDmap)
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
        Nmap(5).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
        
        % calculate a modified MTR map according to RD/BD
        MTR = 100*(PDw-MTw)./(eps+PDw).*(MTw./(eps+PDw)<1.3&MTw./(eps+PDw)>0&PDw>25);
        Nmap(6).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
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
        
        if RFC.RFCorr
            % MFC: We do have P2_a and P2_b parameters for this sequence
            % => T1 = A(B1) + B(B1)*T1app (see Preibisch 2009)
             T1 = RFC.P2_a(1)*f_T.^2+RFC.P2_a(2)*f_T+RFC.P2_a(3)+(RFC.P2_b(1)*f_T.^2+RFC.P2_b(2)*f_T+RFC.P2_b(3)).*((((PDw / fa_pdw) - (T1w / fa_t1w)+eps) ./ ...
                max((T1w * fa_t1w / 2 / TR_t1w) - (PDw * fa_pdw / 2 / TR_pdw),eps))./f_T.^2);
        else
            % MFC: We do not have P2_a or P2_b parameters for this sequence
            % => T1 = T1app
            T1 = ((((PDw / fa_pdw) - (T1w / fa_t1w)+eps) ./ max((T1w * fa_t1w / 2 / TR_t1w) - (PDw * fa_pdw / 2 / TR_pdw),eps))./f_T.^2);
        end
        
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
    Nmap(1).dat(:,:,p) = min(max(tmp,-threshall.R1),threshall.R1); %Truncating images to dynamic range. Originally 20000. Negative values are allowed or min(max(tmp,0),2000)
    
    % A values proportional to PD
    if (~isempty(f_T)) && (~isempty(f_R))
        % Transmit and receive bias corrected quantitative A values
        % again: correct A for transmit bias f_T and receive bias f_R
        % Acorr = A / f_T / f_R , proportional PD
        A = (T1 .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w))./f_T./f_R;
    elseif(~isempty(f_T))&&(isempty(f_R))&&(PDproc.PDmap)
        A = T1 .* (T1w .*(fa_t1w*f_T) / 2 / TR_t1w) + (T1w ./ (fa_t1w*f_T));
    else
        % semi-quantitative A
        A = T1 .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w);
    end
    
    A_forMT = T1_forMT .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w);
    tmp      = A;
    Nmap(2).dat(:,:,p) = max(min(tmp,threshall.A),-threshall.A); % dynamic range increased to 10^5 to accommodate phased-array coils and symmetrical for noise distribution
    
    % MT in [p.u.]; offset by - famt * famt / 2 * 100 where MT_w = 0 (outside mask)
    MT       = ( (A_forMT * fa_mtw - MTw) ./ (MTw+eps) ./ (T1_forMT + eps) * TR_mtw - fa_mtw * fa_mtw / 2 ) * 100;
    if (~isempty(f_T))
        MT = MT .* (1 - 0.4) ./ (1 - 0.4 * f_T);
    end
    if (PDproc.PDmap)%The two outer voxels in all directions are nulled in order to remove artefactual effects from the MT map
        MTforA=MT;
        % Set outer five planes to 0 to prevent errors with bias field
        % estimation during segmentation:
        if (p < 5) || (p > dm(3) - 5)
            MTforA=zeros(size(MT,1),size(MT,2));
        else
            MTforA(1:5,:)=0; MTforA(end-5:end, :)=0;
            MTforA(:, 1:5)=0; MTforA(:, end-5:end)=0;
        end
        tmp      = MTforA;
        if (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw)
            Nmap(7).dat(:,:,p) = max(min(tmp,threshall.MT),-threshall.MT);
        else
            Nmap(5).dat(:,:,p) = max(min(tmp,threshall.MT),-threshall.MT);
        end
    end
    
    tmp      = MT;
    Nmap(3).dat(:,:,p) = max(min(tmp,threshall.MT),-threshall.MT);
    
    % calculate synthetic reference signal at trmt and famt using the
    % rational approximation of the Ernst equation
    S_ref      = A_forMT .* fa_mtw * TR_mtw ./ (T1_forMT+eps) ./ ( TR_mtw ./ (T1_forMT+eps) +  fa_mtw * fa_mtw / 2 );
    % MTR_synt = (S_ref ./ MTw - 1) * 100;
    MTR_synt   = (S_ref-MTw) ./ (S_ref+eps) * 100;
    tmp      = MTR_synt;
    Nmap(4).dat(:,:,p) = max(min(tmp,threshall.MTR_synt),-threshall.MTR_synt);
    spm_progress_bar('Set',p);
    
        
end
if(~isempty(f_T))&&(isempty(f_R))&&(PDproc.PDmap)
    PDcalculation(pth)
end
spm_progress_bar('Clear');

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [] = coreg_mt(P_ref, P_src)
% coregisters the structural images
% for MT protocol

for src_nr=1:size(P_src, 1)
    P_src(src_nr,:);
    VG = spm_vol(P_ref);
    VF = spm_vol(P_src(src_nr,:));
    %coregflags.sep = [2 1];
    coregflags.sep = [4 2];
    x = spm_coreg(VG,VF, coregflags);
    %x  = spm_coreg(mireg(i).VG, mireg(i).VF,flags.estimate);
    M  = inv(spm_matrix(x));
    MM = spm_get_space(deblank(VF.fname));
    spm_get_space(deblank(deblank(VF.fname)), M*MM);
end
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [] = coreg_bias_map(P_ref, P_src)
% coregisters the B1 or receive maps with
% the structurals in the MT protocol

P_src(1,:);
VG = spm_vol(P_ref);
VF = spm_vol(P_src(1,:));
%coregflags.sep = [2 1];
coregflags.sep = [4 2];
x = spm_coreg(VG,VF, coregflags);
%x  = spm_coreg(mireg(i).VG, mireg(i).VF,flags.estimate);
M  = inv(spm_matrix(x));
MM = spm_get_space(deblank(VF.fname));
spm_get_space(deblank(deblank(VF.fname)), M*MM);

VF2 = spm_vol(P_src(2,:)); % now also apply transform to the map
M  = inv(spm_matrix(x));
MM = spm_get_space(deblank(VF2.fname));
spm_get_space(deblank(deblank(VF2.fname)), M*MM);

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function PDcalculation(pth)
disp('----- Calculating Proton Density map -----');

% get PD processing default settings
PDproc = vbq_get_defaults('PDproc');
threshA = vbq_get_defaults('qMRI_maps_thresh.A');

% Creation of whole-brain and white-matter masks
P=spm_select('FPList',pth ,'^.*_MTforA.(img|nii)$');
clear matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel.vols = {P};
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0];
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

TPMs=spm_read_vols(spm_vol(spm_select('FPList',fileparts(P),'^c.*\.(img|nii)$')));
WBmask=zeros(size(squeeze(TPMs(:,:,:,1))));
WBmask(sum(cat(4,TPMs(:,:,:,1:2),TPMs(:,:,:,end)),4)>=PDproc.WBMaskTh)=1;
WMmask=zeros(size(squeeze(TPMs(:,:,:,1))));
WMmask(squeeze(TPMs(:,:,:,2))>=PDproc.WMMaskTh)=1;

temp=spm_select('FPList',pth,'^.*_MTforA');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end
% End of creation of whole-brain and white-matter masks


% Saves masked A map for bias-field correction later
P=spm_select('FPList',pth ,'^.*_A.(img|nii)$');
Vsave=spm_vol(P);
Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),['masked_' spm_str_manip(Vsave.fname,'t')]);
Amap=spm_read_vols(spm_vol(P)).*WBmask;
Amap(Amap==Inf)=0;Amap(isnan(Amap))=0;Amap(Amap==threshA)=0;
spm_write_vol(Vsave,Amap);

% Bias-field correction of masked A map
P=spm_select('FPList',pth ,'^masked.*_A.(img|nii)$');
clear matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel.vols = {P};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = PDproc.biasreg;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = PDproc.biasfwhm;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0];
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

temp=spm_select('FPList',pth,'^(c|masked).*\_A');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end

% Bias field correction of A map. The bias field is calculated on
% the masked A map but we want to apply it on the unmasked A map. We
% therefore need to explicitly load the bias field and apply it on the original A map instead of just
% loading the bias-field corrected A map from the previous step
P=spm_select('FPList',pth ,'^s.*_A.(img|nii)$');
bf = fullfile(pth, spm_select('List', pth, '^BiasField.*\.(img|nii)$'));
BF = double(spm_read_vols(spm_vol(bf)));
Y = BF.*spm_read_vols(spm_vol(P));

% Calibration of flattened A map to % water content using typical white
% matter value from the litterature (69%)

A_WM=WMmask.*Y;
Y=Y/mean(A_WM(A_WM~=0))*69;
sprintf('mean White Matter intensity: %04d',mean(A_WM(A_WM~=0)))
sprintf('SD White Matter intensity %04d',std(A_WM(A_WM~=0),[],1))
Y(Y>200)=0;
% MFC: Estimating Error for data set to catch bias field issues:
errorEstimate = std(A_WM(A_WM > 0))./mean(A_WM(A_WM > 0));
Vsave=spm_vol(P);
Vsave.descrip = ['A Map.  Error Estimate: ', num2str(errorEstimate)];
if errorEstimate > 0.06
    % MFC: Testing on 15 subjects showed 6% is a good cut-off:
    warning(['Error estimate is high: ', Vsave.fname]);
end

% V.fname = P;
spm_write_vol(Vsave,Y);

temp=spm_select('FPList',pth,'^Bias.*\_A');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function bl = feq(val, comp_val)
% floating point comparison
bl = abs(val - comp_val) <= eps(comp_val);

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function p = hinfo(P)
N = nifti(P);
p(numel(N)) = struct('tr',[],'te',[],'fa',[]);
for ii = 1:numel(N),
    tmp = regexp(N(ii).descrip,...
        'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg',...
        'names');
    p(ii).tr = str2num(tmp.tr); %#ok<*ST2NM>
    p(ii).te = str2num(tmp.te);
    p(ii).fa = str2num(tmp.fa);
end
end
