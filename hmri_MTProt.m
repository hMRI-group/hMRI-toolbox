function [fR1, fR2s, fMT, fA, PPDw, PT1w]  = hmri_MTProt(P_mtw, P_pdw, P_t1w, P_trans, P_receiv)

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
%
% MFC 23.10.2015    Adding OLS R2* map option. For details see Weiskopf et 
%                   al., Front. Neurosci. 2014 DOI: 10.3389/fnins.2014.00278
%                   This reference should be cited if you use this output.

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
p = hmri_hinfo(P_mtw);
if ~isempty(p)
    hmri_get_defaults('MPMacq.TE_mtw', cat(1,p.te));
    hmri_get_defaults('MPMacq.TR_mtw', p(1).tr);
    hmri_get_defaults('MPMacq.fa_mtw', p(1).fa);
else
    warning('No TE, TR, and FA values found for MTw images. Fallback to defaults.')
end

p = hmri_hinfo(P_pdw);
if ~isempty(p)
    hmri_get_defaults('MPMacq.TE_pdw', cat(1,p.te));
    hmri_get_defaults('MPMacq.TR_pdw', p(1).tr);
    hmri_get_defaults('MPMacq.fa_pdw', p(1).fa);
else
    warning('No TE, TR, and FA values found for PDw images. Fallback to defaults.')
end

p = hmri_hinfo(P_t1w);
if ~isempty(p)
    hmri_get_defaults('MPMacq.TE_t1w', cat(1,p.te));
    hmri_get_defaults('MPMacq.TR_t1w', p(1).tr);
    hmri_get_defaults('MPMacq.fa_t1w', p(1).fa);
else
    warning('No TE, TR, and FA values found for T1w images. Fallback to defaults.')    
end

json = hmri_get_defaults('json');

% retrieve acquisition parameters
MPMacq = hmri_get_defaults('MPMacq');
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
    
%  Get the MPMacq parameters specific for each protocols:  [TR_pdw TR_t1w fa_pdw fa_t1w]
%  and their specific protocol values/names/tags.
MPMacq_prot = [ ...
hmri_get_defaults('MPMacq.TR_pdw') ...
hmri_get_defaults('MPMacq.TR_t1w') ...
hmri_get_defaults('MPMacq.fa_pdw') ...
hmri_get_defaults('MPMacq.fa_t1w')];

MPMacq_sets = hmri_get_defaults('MPMacq_set');
% then match the values and find protocol tag
Nb_protocols = numel(MPMacq_sets.vals);
ii = 1; mtch = false;
while ~mtch && ii <= Nb_protocols
    if all(MPMacq_prot == MPMacq_sets.vals{ii})
        mtch  = true;
    else
        ii = ii+1;
    end
end
% Define the protocol tag
if ~mtch
    prot_tag = 'Unknown';
else
    prot_tag = MPMacq_sets.tags{ii};
end
% Set the tag for the MPMacq set.
hmri_get_defaults('MPMacq.tag',prot_tag);

%% locally retrieves default parameters from hmri_defaults
% load threshold to save qMRI maps
threshall = hmri_get_defaults('qMRI_maps_thresh');
% load PD maps processing parameters
PDproc = hmri_get_defaults('PDproc');

% RF spoiling correction parameters
RFC = hmri_get_defaults(['rfcorr.',prot_tag]);

% a few words to summarize the situation...
disp(['INFO: Acquisition protocol = ' RFC.tag]);

% check that echo times are identical
for nr = 1:min([length(TE_mtw), length(TE_pdw), length(TE_t1w)])
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

tmp = hmri_get_defaults('outdir');
if ~strcmp(pth, tmp)
  pth = tmp;
end
qMRIcalc=hmri_get_defaults('qMRI_maps');
if (qMRIcalc.QA||PDproc.PDmap||qMRIcalc.ACPCrealign)
    MPMcalcFolder=fullfile(pth,'MPMcalcFolder');
    if(exist(MPMcalcFolder,'dir'))%in case a previous run was stopped prematurely
        rmdir(MPMcalcFolder,'s')
    end
    mkdir(MPMcalcFolder)
end


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
for p = 1:dm(3)
    M = spm_matrix([0 0 p 0 0 0 1 1 1]);
    Y = zeros(dm(1:2));
    for i = 1:n
        M1 = V_pdw(i).mat\V_pdw(1).mat*M;
        Y  = Y + W(i)*log(max(spm_slice_vol(V_pdw(i),M1,dm(1:2),1),1));
    end
    Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'map creation';
Output_hdr.history.procstep.procpar = struct('threshold',threshall);
for ctr = 1:numel(V_pdw)
    Output_hdr.history.input{ctr}.filename = V_pdw(ctr).fname;
    input_hdr = get_metadata(V_pdw(ctr).fname);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{ctr}.history = 'No history available.';
    end
%     Output_hdr.history.input{ctr}.history=input_hdr{1}.history;
end
Output_hdr.history.output.imtype = 'R2* map';
Output_hdr.history.output.units = 'ms-1';
set_metadata(fR2s,Output_hdr,json);

% Average first few echoes for increased SNR and fit T2*
disp('----- Reading and averaging the images -----');

nr_TE_limit = find(TE_mtw > TE_limit,1);
avg_nr      = min([nr_c_echoes nr_TE_limit]);
PDproc.nr_echoes_forA
PP   = {P_mtw,P_pdw,P_t1w};
nam1 = {'MTw','PDw','T1w'};
avg  = [0 0 0];
for ii=1:3
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
    for p=1:dm(3)
        M = spm_matrix([0 0 p]);
        Y = zeros(dm(1:2));
        for nr=1:avg_nr
            M1 = V(nr).mat\V(1).mat*M;
            Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),1);
        end
        Ni.dat(:,:,p) = Y/avg_nr;
        sm = sm + sum(Y(:))/avg_nr;
        spm_progress_bar('Set',p);
    end
    avg(ii) = sm/prod(dm);
    spm_progress_bar('Clear');
    
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = Ni.descrip;
    Output_hdr.history.procstep.procpar = struct('threshold',threshall);
    for ctr = 1:numel(V)
        Output_hdr.history.input{ctr}.filename = V(ctr).fname;
        input_hdr = get_metadata(V(ctr).fname);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
    end
    Output_hdr.history.output.imtype = Ni.descrip;
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(fullfile(pth,[nam '_' nam1{ii} '.nii']),Output_hdr,json);
end
V         = spm_vol(PP{3});
dm        = V(1).dim;
Ni        = nifti;
Ni.mat    = V(1).mat;
Ni.mat0   = V(1).mat;
Ni.descrip= sprintf('Averaged %s images', 'T1w_forA');
Ni.dat    = file_array(fullfile(MPMcalcFolder,[nam '_' 'T1w_forA.nii']),dm,dt, 0,1,0);
create(Ni);
spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
sm = 0;
for p=1:dm(3),
    M = spm_matrix([0 0 p]);
    Y = zeros(dm(1:2));
    for nr=1:PDproc.nr_echoes_forA,
        M1 = V(nr).mat\V(1).mat*M;
        Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),1);
    end
    Ni.dat(:,:,p) = Y/PDproc.nr_echoes_forA;
    sm = sm + sum(Y(:))/PDproc.nr_echoes_forA;
    spm_progress_bar('Set',p);
end


fPD = fullfile(pth,[nam '_' nam1{ii} '.nii']); % TL: could this be deleted?

PPDw = fullfile(pth,[nam '_PDw' '.nii']);
PMTw = fullfile(pth,[nam '_MTw' '.nii']);
PT1w = fullfile(pth,[nam '_T1w' '.nii']);
PT1w_forA = fullfile(MPMcalcFolder,[nam '_T1w_forA' '.nii']);

if true
    disp('----- Coregistering the images -----');
    coreg_mt(PPDw, PMTw);
    coreg_mt(PPDw, PT1w);
    coreg_mt(PPDw, PT1w_forA);
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
VT1w_forA = spm_vol(PT1w_forA);

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

%  --- BEING OLS R2s CODE MFC  ---
if hmri_get_defaults('R2sOLS')
    
    % Calculate OLS R2* map from all echoes (ESTATICS, Weiskopf et al. 2014)
    disp('----- Calculation of OLS R2* map -----');
    
    dt        = [spm_type('float32'),spm_platform('bigend')];
    Ni        = nifti;
    Ni.mat    = V.mat;
    Ni.mat0   = V.mat;
    Ni.descrip='OLS R2* map [1/ms]';
    Ni.dat    = file_array(fullfile(pth,[nam '_R2s_OLS' '.nii']),dm,dt, 0,1,0);
    create(Ni);
    
   % Combine the data and echo times:
    TE = [TE_pdw; TE_mtw; TE_t1w];
    
    nPD = numel(TE_pdw);
    nMT = numel(TE_mtw);
    nT1 = numel(TE_t1w);
    nEchoes = nPD + nMT + nT1;
    
    V_contrasts = spm_vol(P_pdw);
    V_contrasts(nPD+1:nPD+nMT) = spm_vol(P_mtw);
    V_contrasts(nPD+nMT+1:nEchoes) = spm_vol(P_t1w);
    
    % The assumption is that the result of co-registering the average 
    % weighted volumes is applicable for each of the echoes of that
    % contrast => Replicate the mat field across contrasts for all echoes.
    matField = cat(3, repmat(VPDw.mat, [1, 1, nPD]), ...
        repmat(VMTw.mat, [1, 1, nMT]), repmat(VT1w.mat, [1, 1, nT1]));
    
    % Same formalism as for PDw fit but now extra colums for the "S(0)" 
    % amplitudes of the different contrasts:
    reg = [zeros(nEchoes,3) TE(:)];
    reg(1 : nPD, 1)             = 1;
    reg(nPD + 1 : nPD+nMT, 2)   = 1;
    reg(nPD+nMT+1:nEchoes, 3)   = 1;
    W   = (reg'*reg)\reg';
    
    spm_progress_bar('Init',dm(3),'OLS R2* fit','planes completed');
    for p = 1:dm(3),
        M = spm_matrix([0 0 p 0 0 0 1 1 1]);
        data = zeros([nEchoes dm(1:2)]);
        
        for e = 1:nEchoes
            % Take slice p (defined in M) and map to a location in the 
            % appropriate contrast using the matField entry for that
            % contrast, which has been co-registered to the PD-weighted 
            % data:
            M1 = matField(:,:,e)\V_contrasts(1).mat*M;

            % Third order B-spline interpolation for OLS R2* estimation
            % since we no longer assume that the echoes are perfectly 
            % aligned as we do for the standard PDw derived R2* estimate.
            data(e,:,:) = log(max(spm_slice_vol(V_contrasts(e),M1,dm(1:2),3),eps));
        end
        Y = W*reshape(data, [nEchoes prod(dm(1:2))]);
        Y = -reshape(Y(4,:), dm(1:2));
        
        % Written out in Ni with mat field defined by V_templ => first PDw
        % echo:
        Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');
    Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
    Output_hdr.history.procstep.version = hmri_get_version;
    Output_hdr.history.procstep.descrip = 'map creation';
    Output_hdr.history.procstep.procpar = struct('threshold',threshall);
    for ctr = 1:numel(V_pdw)
        Output_hdr.history.input{ctr}.filename = V_pdw(ctr).fname;
        input_hdr = get_metadata(V_pdw(ctr).fname);
        if ~isempty(input_hdr{1})
            Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
        else
            Output_hdr.history.input{ctr}.history = 'No history available.';
        end
        %         Output_hdr.history.input{ctr}.history=input_hdr{1}.history;
    end
    Output_hdr.history.output.imtype = 'R2*-OLS map';
    Output_hdr.history.output.units = 'ms-1';
    set_metadata(fullfile(pth,[nam '_R2s_OLS' '.nii']),Output_hdr,json);
        
end % OLS code

nam2    = {'R1','A','MT','MTR_synt'};
descrip = {'R1 map [1000/s]', 'A map','Delta MT map', 'Synthetic MTR image'};
units = {'1000/s', '%','A.U.', 'A.U.'};
if (TR_mtw == TR_pdw) & (fa_mtw == fa_pdw),
    nam2    = {nam2{:}, 'MTR','MTRdiff'};
    descrip = {descrip{:}, 'Classic MTR image','Percent diff. MTR image (RD/BD)'};
    units = {units{:}, 'A.U.','A.U.'};
end

Nmap    = nifti;
for ii=1:numel(nam2)
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

for p = 1:dm(3)
    M = M0*spm_matrix([0 0 p]);

    MTw = spm_slice_vol(VMTw,VMTw.mat\M,dm(1:2),3);
    PDw = spm_slice_vol(VPDw,VPDw.mat\M,dm(1:2),3);
    T1w = spm_slice_vol(VT1w,VT1w.mat\M,dm(1:2),3);
    T1w_forA = spm_slice_vol(VT1w_forA,VT1w_forA.mat\M,dm(1:2),3);
    
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
             T1 = RFC.P2_a(1)*f_T.^2 + ...
                  RFC.P2_a(2)*f_T + ...
                  RFC.P2_a(3) + ...
                  (RFC.P2_b(1)*f_T.^2+RFC.P2_b(2)*f_T+RFC.P2_b(3)) .* ...
                  ((((PDw / fa_pdw) - (T1w / fa_t1w)+eps) ./ ...
                  max((T1w * fa_t1w / 2 / TR_t1w) - (PDw * fa_pdw / 2 / TR_pdw),eps))./f_T.^2);
        else
            % MFC: We do not have P2_a or P2_b parameters for this sequence
            % => T1 = T1app
            T1 = ((((PDw / fa_pdw) - (T1w / fa_t1w)+eps) ./ ...
                max((T1w * fa_t1w / 2 / TR_t1w) - (PDw * fa_pdw / 2 / TR_pdw),eps))./f_T.^2);
        end
        
        R1 = 1./T1*10^6;
        
        %          R1App=f_T.^2.*(((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw)) ./ ...
        %              max(((PDw / fa_pdw) - (T1w / fa_t1w)),eps));
        %          A_poly=P2_a(1)*f_T.^2+P2_a(2)*f_T+P2_a(3);B_poly=(P2_b(1)*f_T.^2+P2_b(2)*f_T+P2_b(3));
        %          R1 = R1App./(max(R1App.*A_poly+B_poly,eps))*10^6;
    end
    T1_forMT = ((PDw / fa_pdw) - (T1w / fa_t1w)) ./ ...
        max((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw),eps);
    T1       = max(T1,0);

    
    R1(R1<0) = 0;
    tmp      = R1;
    Nmap(1).dat(:,:,p) = min(max(tmp,-threshall.R1),threshall.R1);
    % Truncating images to dynamic range. Originally 20000. Negative values are allowed or min(max(tmp,0),2000)
    
    % A values proportional to PD
    if (~isempty(f_T)) && (~isempty(f_R))
        % Transmit and receive bias corrected quantitative A values
        % again: correct A for transmit bias f_T and receive bias f_R
        % Acorr = A / f_T / f_R , proportional PD
        A = (T1 .* (T1w_forA * fa_t1w / 2 / TR_t1w) + (T1w_forA / fa_t1w))./f_T./f_R;
    elseif(~isempty(f_T))&&(isempty(f_R))%&&(PDproc.PDmap)
        A = T1 .* (T1w_forA .*(fa_t1w*f_T) / 2 / TR_t1w) + (T1w_forA ./ (fa_t1w*f_T));
    else
        % semi-quantitative A
        A = T1 .* (T1w_forA * fa_t1w / 2 / TR_t1w) + (T1w_forA / fa_t1w);
    end
    
    A_forMT = T1_forMT .* (T1w * fa_t1w / 2 / TR_t1w) + (T1w / fa_t1w);
    tmp      = A;
    Nmap(2).dat(:,:,p) = max(min(tmp,threshall.A),-threshall.A);
    % dynamic range increased to 10^5 to accommodate phased-array coils and symmetrical for noise distribution
    
    % MT in [p.u.]; offset by - famt * famt / 2 * 100 where MT_w = 0 (outside mask)
    MT       = ( (A_forMT * fa_mtw - MTw) ./ (MTw+eps) ./ (T1_forMT + eps) * TR_mtw - fa_mtw * fa_mtw / 2 ) * 100;
    if (~isempty(f_T))
        MT = MT .* (1 - 0.4) ./ (1 - 0.4 * f_T);
    end
    
    tmp      = MT;
    Nmap(3).dat(:,:,p) = max(min(tmp,threshall.MT),-threshall.MT);
    
    % calculate synthetic reference signal at trmt and famt using the
    % rational approximation of the Ernst equation
    S_ref      = A_forMT .* fa_mtw * TR_mtw ./ (T1_forMT+eps) ./ ...
                ( TR_mtw ./ (T1_forMT+eps) +  fa_mtw * fa_mtw / 2 );
    % MTR_synt = (S_ref ./ MTw - 1) * 100;
    MTR_synt   = (S_ref-MTw) ./ (S_ref+eps) * 100;
    tmp      = MTR_synt;
    Nmap(4).dat(:,:,p) = max(min(tmp,threshall.MTR_synt),-threshall.MTR_synt);
    spm_progress_bar('Set',p);
end

Vtemp = cat(1,VMTw,VPDw,VT1w);
Output_hdr = struct('history',struct('procstep',[],'input',[],'output',[]));
Output_hdr.history.procstep.version = hmri_get_version;
Output_hdr.history.procstep.descrip = 'map creation';
Output_hdr.history.procstep.procpar = struct('threshold',threshall);
for ctr = 1:numel(Vtemp)
    Output_hdr.history.input{ctr}.filename = Vtemp(ctr).fname;
    input_hdr = get_metadata(Vtemp(ctr).fname);
    if ~isempty(input_hdr{1})
        Output_hdr.history.input{ctr}.history = input_hdr{1}.history;
    else
        Output_hdr.history.input{ctr}.history = 'No history available.';
    end
end
for ctr = 1:size(nam2,2)
    Output_hdr.history.output.imtype = descrip(ctr);
    Output_hdr.history.output.units = units(ctr);
    set_metadata(fullfile(pth,[nam '_' nam2{ctr} '.nii']),Output_hdr,json);
end

if (qMRIcalc.QA||(PDproc.PDmap))
    R = spm_select('FPList',pth,'^s.*_MT.(img|nii)$');
    copyfile(R(1,:),MPMcalcFolder);
    Vsave = spm_vol(spm_select('FPList',MPMcalcFolder,'^s.*_MT.(img|nii)$'));
    MTtemp=spm_read_vols(Vsave);
    %The 5 outer voxels in all directions are nulled in order to remove artefactual effects from the MT map on segmentation:
    MTtemp(1:5,:,:)=0; MTtemp(end-5:end, :,:)=0;
    MTtemp(:, 1:5,:)=0; MTtemp(:, end-5:end,:)=0;
    MTtemp(:,:, 1:5)=0; MTtemp(:,:,end-5:end)=0;
    spm_write_vol(Vsave,MTtemp);
    
    P = spm_select('FPList',MPMcalcFolder,'^s.*_MT.(img|nii)$');
    clear matlabbatch
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {P};
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
%     spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);
end

if ~isempty(f_T) && isempty(f_R) && PDproc.PDmap
%     PDcalculation(pth)
    PDcalculation(pth,MPMcalcFolder)
end
if (qMRIcalc.QA||PDproc.PDmap||qMRIcalc.ACPCrealign)
    rmdir(MPMcalcFolder,'s')
end

spm_progress_bar('Clear');

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [] = coreg_mt(P_ref, P_src)
% coregisters the structural images
% for MT protocol

for src_nr = 1:size(P_src,1)
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
function PDcalculation(pth,MPMcalcFolder)
disp('----- Calculating Proton Density map -----');

PDproc = hmri_get_defaults('PDproc');
threshA = hmri_get_defaults('qMRI_maps_thresh.A');

TPMs=spm_read_vols(spm_vol(spm_select('FPList',MPMcalcFolder,'^c.*\.(img|nii)$')));
WBmask=zeros(size(squeeze(TPMs(:,:,:,1))));
WBmask(sum(cat(4,TPMs(:,:,:,1:2),TPMs(:,:,:,end)),4)>=PDproc.WBMaskTh)=1;
WMmask=zeros(size(squeeze(TPMs(:,:,:,1))));
WMmask(squeeze(TPMs(:,:,:,2))>=PDproc.WBMaskTh)=1;

% Saves masked A map for bias-field correction later
P=spm_select('FPList',pth ,'^.*_A.(img|nii)$');
Vsave=spm_vol(P);
Vsave.fname=fullfile(MPMcalcFolder,['masked_' spm_str_manip(Vsave.fname,'t')]);
Amap=spm_read_vols(spm_vol(P)).*WBmask;
Amap(Amap==Inf)=0;Amap(isnan(Amap))=0;Amap(Amap==threshA)=0;
spm_write_vol(Vsave,Amap);

% Bias-field correction of masked A map
P=spm_select('FPList',MPMcalcFolder ,'^masked.*_A.(img|nii)$');
clear matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel.vols = {P};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = PDproc.biasreg;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = PDproc.biasfwhm;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0];
% spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

temp=spm_select('FPList',MPMcalcFolder,'^(c|masked).*\_A');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end

% Bias field correction of A map. The bias field is calculated on
% the masked A map but we want to apply it on the unmasked A map. We
% therefore need to explicitly load the bias field and apply it on the original A map instead of just
% loading the bias-field corrected A map from the previous step
P=spm_select('FPList',pth ,'^s.*_A.(img|nii)$');
bf = fullfile(MPMcalcFolder, spm_select('List', MPMcalcFolder, '^BiasField.*\.(img|nii)$'));
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
if vbq_get_defaults('QA')
    if (exist(fullfile(pth,'QualityAssessment.mat'),'file')==2)
        load(fullfile(pth,'QualityAssessment.mat'));
    end
    QA.PD.mean=mean(A_WM(A_WM > 0));QA.PD.SD=std(A_WM(A_WM > 0));
    save(fullfile(pth,'QualityAssessment'),'QA')
end
% V.fname = P;
spm_write_vol(Vsave,Y);

temp=spm_select('FPList',MPMcalcFolder,'^Bias.*\_A');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function bl = feq(val, comp_val)
% floating point comparison
bl = abs(val - comp_val) <= eps(comp_val);

end
