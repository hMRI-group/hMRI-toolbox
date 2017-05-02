function [fR1, fR2s, fMT, fA, PPDw, PT1w]  = hmri_MTProt(jobsubj, P_trans, P_receiv) %#ok<*STOUT>

% Evaluation function for multi-contrast multi-echo FLASH protocol
% P_mtw, P_pdw, P_t1w (retrieved from jobsubj.raw_mpm): MTw, PDw, T1w
%           images (can be multiple echoes = images) 
% TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w: echo times and TR of
%           images 
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


%% =======================================================================%
% Initiate processing
%=========================================================================%
fprintf(1, ['---------------------- MPM PROCESSING ----------------------\n' ...
    'Create maps from multi-contrast multi-echo FLASH protocol.\n']); 
                    
% retrieves all required parameters for MPM processing
mpm_params = get_mpm_params(jobsubj);

% for convenience, define a few parameters to make formulae more readable
% and avoid number of repetitions:
TE_pdw = mpm_params.input.PDw.TE;
TE_mtw = mpm_params.input.MTw.TE;
TE_t1w = mpm_params.input.T1w.TE;
TR_pdw = mpm_params.input.PDw.TR;
TR_mtw = mpm_params.input.MTw.TR;
TR_t1w = mpm_params.input.T1w.TR;
fa_pdw = mpm_params.input.PDw.fa;
fa_mtw = mpm_params.input.MTw.fa;
fa_t1w = mpm_params.input.T1w.fa;
threshall = mpm_params.proc.threshall;
PDproc = mpm_params.proc.PD;
RFC = mpm_params.proc.RFC;
dt = [spm_type('float32'),spm_platform('bigend')]; % for nifti output
outbasename = spm_file(mpm_params.input.MTw.fname(1,:),'basename'); % for all output files
calcpath = mpm_params.calcpath;
mpm_params.outbasename = outbasename;
respath = mpm_params.respath;

% Load B1 mapping data if available 
% P_trans(1,:) = magnitude image (anatomical reference for coregistration) 
% P_trans(2,:) = B1 map (p.u.)
V_trans = [];
if ~isempty(P_trans); V_trans = spm_vol(P_trans); end

% Load sensitivity map if available
% P_receiv(1,:) = magnitude image (anatomical reference for coregistration) 
% P_receiv(2,:) = sensitivity map
V_receiv   = [];
if ~isempty(P_receiv); V_receiv = spm_vol(P_receiv); end


%% =======================================================================%
% Calculate R2* map from PDw echoes
%=========================================================================%
fprintf(1,'\n    -------- R2* map calculation --------\n');

% load PDw images
V_pdw = spm_vol(mpm_params.input.PDw.fname);
dm = V_pdw(1).dim;
spm_progress_bar('Init',dm(3),'R2* fit','planes completed');

% create nifti object for output R2* map
Ni          = nifti;
Ni.mat      = V_pdw(1).mat;
Ni.mat0     = V_pdw(1).mat;
Ni.descrip  = 'R2* map [1/ms]';
Ni.dat      = file_array(fullfile(calcpath,[outbasename '_R2s' '.nii']),dm,dt, 0,1,0);
create(Ni);
fR2s = fullfile(calcpath,[outbasename '_R2s' '.nii']);

% Fit R2* decay
reg = [ones(numel(TE_pdw),1) TE_pdw(:)];
W   = (reg'*reg)\reg';
W   = -W(2,:)';
for p = 1:dm(3)
    M = spm_matrix([0 0 p 0 0 0 1 1 1]);
    Y = zeros(dm(1:2));
    for i = 1:numel(V_pdw)
        M1 = V_pdw(i).mat\V_pdw(1).mat*M;
        Y  = Y + W(i)*log(max(spm_slice_vol(V_pdw(i),M1,dm(1:2),1),1));
    end
    Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

% Set and write metadata
input_files = mpm_params.input.PDw.fname;
Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
Output_hdr.history.output.imtype = 'R2* map';
Output_hdr.history.output.units = 'ms-1';
set_metadata(fR2s,Output_hdr,mpm_params.json);


%% =======================================================================%
% Reading and averaging the images 
%=========================================================================%
fprintf(1,'\n    -------- Reading and averaging the images --------\n');

% Average only first few echoes for increased SNR and fit T2*
nr_TE_limit = sum(TE_mtw < mpm_params.input.TE_limit) + 1; % number of echoes under the TE_limit + 1(not sure why +1, why not to define TE_limit bigger then?)
nr_c_echoes = min([length(TE_mtw), length(TE_pdw), length(TE_t1w)]); % maximum number of echoes available for ALL contrasts
avg_nr      = min([nr_c_echoes nr_TE_limit]); % average is made over maximum number of echoes available for ALL contrasts AND under TE_limit
PP   = {mpm_params.input.MTw.fname,mpm_params.input.PDw.fname,mpm_params.input.T1w.fname}; % gather all images in cell array
contrastnam = {'MTw','PDw','T1w'};
avg  = [0 0 0]; % not used?
for ii=1:3 % loop over MTw, PDw, T1w contrasts
    avg_fnam    = fullfile(calcpath,[outbasename '_' contrastnam{ii} '.nii']);
    eval(sprintf('P%s = avg_fnam;', contrastnam{ii})); % i.e. PPDw/PMTw/PT1w = avg_fnam; Defined here!!
    V           = spm_vol(PP{ii});
    dm          = V(1).dim;
    Ni          = nifti;
    Ni.mat      = V(1).mat;
    Ni.mat0     = V(1).mat;
    Ni.descrip  = sprintf('Averaged %s images', contrastnam{ii});
    Ni.dat      = file_array(avg_fnam,dm,dt,0,1,0);
    create(Ni);
    spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
    sm = 0;
    for p = 1:dm(3)
        M = spm_matrix([0 0 p]);
        Y = zeros(dm(1:2));
        for nr = 1:avg_nr
            M1 = V(nr).mat\V(1).mat*M;
            Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),1);
        end
        Ni.dat(:,:,p) = Y/avg_nr;
        sm = sm + sum(Y(:))/avg_nr;
        spm_progress_bar('Set',p);
    end
    avg(ii) = sm/prod(dm);
    spm_progress_bar('Clear');
    
    input_files = mpm_params.input.(contrastnam{ii}).fname;
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = Ni.descrip;
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(avg_fnam,Output_hdr,mpm_params.json);
end

% Average T1w image for PD calculation (actually, by default, no average
% is performed, only first echo is kept):  
PT1w_forA = fullfile(calcpath,[outbasename '_' 'T1w_forA.nii']);
V           = spm_vol(PP{3});
dm          = V(1).dim;
Ni          = nifti;
Ni.mat      = V(1).mat;
Ni.mat0     = V(1).mat;
Ni.descrip  = 'Averaged T1w images for PD calculation';
Ni.dat      = file_array(PT1w_forA,dm,dt, 0,1,0);
create(Ni);
spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
sm = 0;
for p = 1:dm(3),
    M = spm_matrix([0 0 p]);
    Y = zeros(dm(1:2));
    for nr = 1:PDproc.nr_echoes_forA,
        M1 = V(nr).mat\V(1).mat*M;
        Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),1);
    end
    Ni.dat(:,:,p) = Y/PDproc.nr_echoes_forA;
    sm = sm + sum(Y(:))/PDproc.nr_echoes_forA;
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');


%% =======================================================================%
% Coregistering the images
%=========================================================================%
fprintf(1,'\n    -------- Coregistering the images  --------\n');

% NOTE: PPDw, PMTw and PT1w are defined above as evaluated string (~line 133)
x_MT2PD = coreg_mt(PPDw, PMTw);  %#ok<NODEF>
x_T12PD = coreg_mt(PPDw, PT1w); %#ok<NODEF>
coreg_mt(PPDw, PT1w_forA);
if ~isempty(V_trans)
    coreg_bias_map(PPDw, P_trans);
end
if ~isempty(V_receiv)
    coreg_bias_map(PPDw, P_receiv);
end

% parameters saved for quality assessment
if mpm_params.QA.enable
    if exist(mpm_params.QA.fnam,'file')
        mpm_params.QA = spm_jsonread(mpm_params.QA.fnam);
    end
    mpm_params.QA.ContrastCoreg.MT2PD = x_MT2PD;
    mpm_params.QA.ContrastCoreg.T12PD = x_T12PD;
    spm_jsonwrite(mpm_params.QA.fnam, mpm_params.QA, struct('indent','\t'));
end

% load averaged images
VPDw = spm_vol(PPDw);
VMTw = spm_vol(PMTw);
VT1w = spm_vol(PT1w);
VT1w_forA = spm_vol(PT1w_forA);


%% =======================================================================%
% multi-contrast R2* map calculation for quality assessment by AL
% METHOD: R2s map calculated for each contrast separately (PDw, MTw, T1w)
% to evaluate the SD of each R2* map in white matter (i.e. intra-run motion
% measurement).
%=========================================================================%
if mpm_params.QA.enable
    fprintf(1,'\n    -------- multi-contrast R2* map calculation for QA --------\n');
    
    allTEs = {TE_mtw, TE_pdw, TE_t1w};
    contrastnam = {'MTw', 'PDw', 'T1w'};
    V_all = [VMTw VPDw VT1w];
    V_PD = spm_vol(PP{2});
    for ctr = 1:size(PP,2)
        dt        = [spm_type('float32'),spm_platform('bigend')];
        Ni        = nifti;
        Ni.mat    = V.mat;
        Ni.mat0   = V.mat;
        Ni.descrip='OLS R2* map [1/ms]';
        Ni.dat    = file_array(fullfile(calcpath,[outbasename '_R2s_' contrastnam{ctr} '.nii']),dm,dt, 0,1,0);
        create(Ni);
        
        TE = allTEs{ctr};
        V_contrasts = spm_vol(PP{ctr});
        % The assumption is that the result of co-registering the average
        % weighted volumes is applicable for each of the echoes of that
        % contrast => Replicate the mat field across contrasts for all echoes.
        % % matField = cat(3, repmat(VPDw.mat, [1, 1, nPD]), ...
        % % repmat(VMTw.mat, [1, 1, nMT]), repmat(VT1w.mat, [1, 1, nT1]));

        reg = [ones(size(TE)) TE(:)];
        W   = (reg'*reg)\reg';
        
        spm_progress_bar('Init',dm(3),'multi-contrast R2* fit','planes completed');
        for p = 1:dm(3),
            M = spm_matrix([0 0 p 0 0 0 1 1 1]);
            data = zeros([size(TE,1) dm(1:2)]);
            
            for e = 1:size(TE,1)
                % Take slice p (defined in M) and map to a location in the
                % appropriate contrast using the matField entry for that
                % contrast, which has been co-registered to the PD-weighted
                % data:
                M1 = V_all(ctr).mat\V_PD(1).mat*M;
                
                % Third order B-spline interpolation for OLS R2* estimation
                % since we no longer assume that the echoes are perfectly
                % aligned as we do for the standard PDw derived R2* estimate.
                data(e,:,:) = log(max(spm_slice_vol(V_contrasts(e),M1,dm(1:2),3),eps));
            end
            Y = W*reshape(data, [size(TE,1) prod(dm(1:2))]);
            Y = -reshape(Y(2,:), dm(1:2));
            
            % NB: mat field defined by V_pdw => first PDw echo
            % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec),
            % negative values are allowed to preserve Gaussian distribution
            Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); 
            spm_progress_bar('Set',p);
        end
        spm_progress_bar('Clear');
    end
end

%% =======================================================================%
% OLS R2* map calculation by MFC
% [Reference: ESTATICS, Weiskopf et al. 2014]
%=========================================================================%
if mpm_params.proc.R2sOLS
    fprintf(1,'\n    -------- OLS R2* map calculation --------\n');
    
    R2sOLS_fnam = fullfile(calcpath,[outbasename '_R2s_OLS' '.nii']);
    Ni          = nifti;
    Ni.mat      = V_pdw(1).mat;
    Ni.mat0     = V_pdw(1).mat;
    Ni.descrip  = 'OLS R2* map [1/ms]';
    Ni.dat      = file_array(R2sOLS_fnam,dm,dt,0,1,0);
    create(Ni);
    
    % Combine the data and echo times:
    TE = [TE_pdw; TE_mtw; TE_t1w]; 
    nPD = numel(TE_pdw);
    nMT = numel(TE_mtw);
    nT1 = numel(TE_t1w);
    nEchoes = nPD + nMT + nT1;
    
    % list all volumes PDw, MTw, T1w
    V_contrasts = V_pdw;
    V_contrasts(nPD+1:nPD+nMT) = spm_vol(mpm_params.input.MTw.fname);
    V_contrasts(nPD+nMT+1:nEchoes) = spm_vol(mpm_params.input.T1w.fname);
    
    % The assumption is that the result of co-registering the average 
    % weighted volumes is applicable for each of the echoes of that
    % contrast => Replicate the mat field across contrasts for all echoes.
    matField = cat(3, repmat(VPDw.mat, [1, 1, nPD]), ...
        repmat(VMTw.mat, [1, 1, nMT]), repmat(VT1w.mat, [1, 1, nT1]));
    
    % Same formalism as for PDw fit but now extra colums for the "S(0)" 
    % amplitudes of the different contrasts:
    reg = [zeros(nEchoes,3) TE(:)];
    reg(1 : nPD, 1) = 1;
    reg(nPD + 1 : nPD+nMT, 2) = 1;
    reg(nPD+nMT+1:nEchoes, 3) = 1;
    W = (reg'*reg)\reg';
    
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
        
        % NB: mat field defined by V_pdw => first PDw echo
        % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), 
        % negative values are allowed to preserve Gaussian distribution.
        Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); 
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');

    input_files = mpm_params.input.PDw.fname;
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = 'R2*-OLS map';
    Output_hdr.history.output.units = 'ms-1';
    set_metadata(fullfile(calcpath,[outbasename '_R2s_OLS' '.nii']),Output_hdr,mpm_params.json);
        
end % OLS code

%% =======================================================================%
% Prepare output for R1, PD and MT maps
%=========================================================================%
% description fields and file names of output images
if ~isempty(V_trans) && isempty(V_receiv) && PDproc.PDmap % whether quantitative PD or not...
    output_suffix    = {'R1','PD','MT'};
    descrip = {'R1 map [1000/s]', 'Water concentration [%]','Delta MT map'};
    units = {'1000/s', '%','a.u.'};
else
    output_suffix    = {'R1','A','MT'};
    descrip = {'R1 map [1000/s]', 'Signal amplitude [a.u.]','Delta MT map'};
    units = {'1000/s', 'a.u.','a.u.'};
end
if (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw) % additional MTR image...
    output_suffix    = [output_suffix{:} {'MTR'}];
    descrip = [descrip{:} {'Classic MTR image'}];
    units = [units {'a.u.'}];
end

% define NIFTI objects for output images
Nmap    = nifti;
for ii=1:numel(output_suffix)
    dm        = V_pdw(1).dim;
    Ni        = nifti;
    Ni.mat    = V_pdw(1).mat;
    Ni.mat0   = V_pdw(1).mat;
    Ni.descrip= descrip{ii};
    Ni.dat    = file_array(fullfile(calcpath,[outbasename '_' output_suffix{ii} '.nii']),dm,dt, 0,1,0);
    create(Ni);
    Nmap(ii) = Ni;
end
fR1 = fullfile(calcpath,[outbasename '_' output_suffix{1} '.nii']);
fA  = fullfile(calcpath,[outbasename '_' output_suffix{2} '.nii']);
fMT = fullfile(calcpath,[outbasename '_' output_suffix{3} '.nii']);


%% =======================================================================%
% Map calculation continued (R1, PD, MT) 
%=========================================================================%
fprintf(1,'\n    -------- Map calculation continued (R1, PD, MT) --------\n');

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
    
    if ~isempty(V_trans)
        f_T = spm_slice_vol(V_trans(2,:),V_trans(2,:).mat\M,dm(1:2),3)/100; % divide by 100, since p.u. maps
    else
        f_T = [];
    end
    if ~isempty(V_receiv) && ~isempty(V_trans)
        f_R = spm_slice_vol(V_receiv(2,:),V_receiv(2,:).mat\M,dm(1:2),3)/100; % divide by 100, since p.u. maps
        f_R = f_R .* f_T; % f_R is only the sensitivity map and not the true receive bias map, therefore needs to be multiplied by transmit bias (B1+ approx. B1- map)
    else
        f_R = [];
    end
    
    % Standard magnetization transfer ratio (MTR) in percent units [p.u.]
    % only if  trpd = trmt and fapd = fmt
    % else calculate "synthetic MTR using A and T1 (see below)
    if numel(Nmap)>3 && (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw),
        MTR = (PDw-MTw)./(PDw+eps) * 100;
        % write MTR image
        Nmap(4).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
        
        % % calculate a modified MTR map according to RD/BD
        % MTR = 100*(PDw-MTw)./(eps+PDw).*(MTw./(eps+PDw)<1.3&MTw./(eps+PDw)>0&PDw>25);
        % Nmap(6).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
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
        
        % R1App = f_T.^2.*(((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw)) ./ ...
        %       max(((PDw / fa_pdw) - (T1w / fa_t1w)),eps));
        % A_poly = P2_a(1)*f_T.^2+P2_a(2)*f_T+P2_a(3);B_poly=(P2_b(1)*f_T.^2+P2_b(2)*f_T+P2_b(3));
        % R1 = R1App./(max(R1App.*A_poly+B_poly,eps))*10^6;
    end
    T1_forMT = ((PDw / fa_pdw) - (T1w / fa_t1w)) ./ ...
        max((T1w * (fa_t1w / 2 / TR_t1w)) - (PDw * fa_pdw / 2 / TR_pdw),eps);
    T1       = max(T1,0);
    
    R1(R1<0) = 0;
    tmp      = R1;
    Nmap(1).dat(:,:,p) = min(max(tmp,-threshall.R1),threshall.R1); % truncating images
    
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
    
    % % calculate synthetic reference signal at trmt and famt using the
    % % rational approximation of the Ernst equation
    % S_ref = A_forMT .* fa_mtw * TR_mtw ./ (T1_forMT+eps) ./ ...
    %       ( TR_mtw ./ (T1_forMT+eps) +  fa_mtw * fa_mtw / 2 );
    % % MTR_synt = (S_ref ./ MTw - 1) * 100;
    % MTR_synt = (S_ref-MTw) ./ (S_ref+eps) * 100;
    % tmp = MTR_synt;
    % Nmap(4).dat(:,:,p) = max(min(tmp,threshall.MTR_synt),-threshall.MTR_synt);
    % spm_progress_bar('Set',p);
end

% set metadata for all output images
input_files = cat(1,mpm_params.input.MTw.fname,mpm_params.input.PDw.fname,mpm_params.input.T1w.fname);
Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
for ctr = 1:size(output_suffix,2)
    Output_hdr.history.output.imtype = descrip(ctr);
    Output_hdr.history.output.units = units(ctr);
    set_metadata(fullfile(calcpath,[outbasename '_' output_suffix{ctr} '.nii']),Output_hdr,mpm_params.json);
end

%% =======================================================================%
% ACPC Realign all images
%=========================================================================%
if mpm_params.ACPCrealign
    fprintf(1,'\n    -------- ACPC Realign all images --------\n');
    
    % Define and calculate masked MT image
    % Load MT image
    V_MT = spm_vol(fMT);
    MTimage = spm_read_vols(V_MT);    
    % Define new file name for masked MT image
    V_MT.fname = fullfile(calcpath,['masked_' spm_str_manip(fMT,'t')]);
    % Load average PDw image (mask based on averaged PDw image)
    PDWimage = spm_read_vols(spm_vol(PPDw));
    % Mask MT image and save the masked MT image ('masked_..._MT.nii')
    MTimage(PDWimage<0.6*mean(PDWimage(:)))=0;
    spm_write_vol(V_MT,MTimage);

    % Use masked MT image to calculate transformation for ACPC realignment
    % (to increase robustness in segmentation):
    [~,R] = hmri_comm_adjust(1,V_MT.fname,V_MT.fname,8,0,fullfile(spm('Dir'),'canonical','avg152T1.nii')); 
    
    % Collect all images from all defined output directories:
    allpaths = jobsubj.path;
    ACPC_images = [];
    f = fieldnames(allpaths);
    for cfield = 1:length(f)
        tmpfiles = spm_select('FPList',allpaths.(f{cfield}),'.*.(img|nii)$');
        if ~isempty(tmpfiles)
            if ~isempty(ACPC_images)
                ACPC_images = char(ACPC_images,spm_select('FPList',allpaths.(f{cfield}),'.*.(img|nii)$'));
            else
                ACPC_images = tmpfiles;
            end
        end
    end
    
    % Apply transform to ALL images 
    for i=1:size(ACPC_images,1)
        spm_get_space(deblank(ACPC_images(i,:)),...
            R*spm_get_space(deblank(ACPC_images(i,:))));
        Vsave = spm_vol(ACPC_images(i,:));
        Vsave.descrip = [Vsave.descrip ' - AC-PC realigned'];
        spm_write_vol(Vsave,spm_read_vols(spm_vol(ACPC_images(i,:))));
    end;
    
    % Save transformation matrix
    spm_jsonwrite(fullfile(respath,'MPM_map_creation_ACPCrealign_transformation_matrix.json'),R,struct('indent','\t'));
end

% for quality assessment and/or PD map calculation
if (mpm_params.QA.enable||(PDproc.PDmap))
    Vsave = spm_vol(fMT);
    MTtemp = spm_read_vols(Vsave);
    %The 5 outer voxels in all directions are nulled in order to remove artefactual effects from the MT map on segmentation:
    MTtemp(1:5,:,:)=0; MTtemp(end-5:end,:,:)=0;
    MTtemp(:,1:5,:)=0; MTtemp(:,end-5:end,:)=0;
    MTtemp(:,:,1:5)=0; MTtemp(:,:,end-5:end)=0;
    spm_write_vol(Vsave,MTtemp);
    
    clear matlabbatch
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {fMT};
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    spm_jobman('run', matlabbatch);
end

% for quality assessment
if mpm_params.QA.enable
    % Calculate WM mask
    TPMs = spm_read_vols(spm_vol(spm_select('FPList',calcpath,'^c.*\.(img|nii)$')));
    WMmask = zeros(size(squeeze(TPMs(:,:,:,2))));
    WMmask(squeeze(TPMs(:,:,:,2))>=PDproc.WMMaskTh) = 1;        
    WMmask = spm_erode(spm_erode(double(WMmask)));
    % Load OLS R2s maps calculated for each contrast and mask tehm
    R2s = spm_read_vols(spm_vol(spm_select('FPList',calcpath,'^s.*_R2s_(MTw|PDw|T1w).(img|nii)$')));
    R2s = R2s.*repmat(WMmask,[1 1 1 size(R2s,4)]);
    SDR2s = zeros(1,size(R2s,4));
    % For each contrast calculate SD of the R2s values within the WM mask
    % (measure of the intra-run motion for each contrast)
    for ctr=1:size(R2s,4)
        MaskedR2s = squeeze(R2s(:,:,:,ctr));
        SDR2s(ctr) = std(MaskedR2s(MaskedR2s~=0),[],1);
    end
    if exist(mpm_params.QA.fnam,'file')
        mpm_params.QA = spm_jsonread(mpm_params.QA.fnam);
    end
    mpm_params.QA.SDR2s.MTw = SDR2s(1);
    mpm_params.QA.SDR2s.PDw = SDR2s(2);
    mpm_params.QA.SDR2s.T1w = SDR2s(3);
    spm_jsonwrite(mpm_params.QA.fnam, mpm_params.QA, struct('indent','\t'));
end

% PD map calculation
if ~isempty(f_T) && isempty(f_R) && PDproc.PDmap
    PDcalculation(fA, mpm_params);
end

% copy final result files into Results directory
fR1_final = fullfile(respath, spm_file(fR1,'filename'));
copyfile(fR1,fR1_final);
try copyfile([spm_str_manip(fR1,'r') '.json'],[spm_str_manip(fR1_final,'r') '.json']); end %#ok<*TRYNC>
fR1 = fR1_final;

fR2s_final = fullfile(respath, spm_file(fR2s,'filename'));
copyfile(fR2s,fR2s_final);
try copyfile([spm_str_manip(fR2s,'r') '.json'],[spm_str_manip(fR2s_final,'r') '.json']); end
fR2s = fR2s_final;

fMT_final = fullfile(respath, spm_file(fMT,'filename'));
copyfile(fMT,fMT_final);
try copyfile([spm_str_manip(fMT,'r') '.json'],[spm_str_manip(fMT_final,'r') '.json']); end
fMT = fMT_final;

fA_final = fullfile(respath, spm_file(fA,'filename'));
copyfile(fA,fA_final);
try copyfile([spm_str_manip(fA,'r') '.json'],[spm_str_manip(fA_final,'r') '.json']); end
fA = fA_final;

PPDw_final = fullfile(respath, spm_file(PPDw,'filename'));
copyfile(PPDw,PPDw_final);
try copyfile([spm_str_manip(PPDw,'r') '.json'],[spm_str_manip(PPDw_final,'r') '.json']); end
PPDw = PPDw_final;

PT1w_final = fullfile(respath, spm_file(PT1w,'filename'));
copyfile(PT1w,PT1w_final);
try copyfile([spm_str_manip(PT1w,'r') '.json'],[spm_str_manip(PT1w_final,'r') '.json']); end
PT1w = PT1w_final;

% save processing params (mpm_params)
spm_jsonwrite(fullfile(respath,'MPM_map_creation_mpm_params.json'),mpm_params,struct('indent','\t'));

spm_progress_bar('Clear');

end

%% =======================================================================%
% To coregister the structural images for MT protocol 
%=========================================================================%
function [x] = coreg_mt(P_ref, P_src)

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

%% =======================================================================%
% To coregister the B1 or receive maps with the anatomical images in the
% MT protocol. 
%=========================================================================%
function [] = coreg_bias_map(P_ref, P_src)

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

%% =======================================================================%
% Proton density map calculation
%=========================================================================%
function PDcalculation(fA, mpm_params)
% fA is the filename of the output A image 
% (not yet properly quantitative PD map when it enters PDcalculation)
fprintf(1,'\n    -------- Proton density map calculation --------\n');

PDproc = mpm_params.proc.PD;
threshA = mpm_params.proc.threshall.A;
calcpath = mpm_params.calcpath;

TPMs = spm_read_vols(spm_vol(spm_select('FPList',calcpath,'^c.*\.(img|nii)$')));
WBmask = zeros(size(squeeze(TPMs(:,:,:,1))));
WBmask(sum(cat(4,TPMs(:,:,:,1:2),TPMs(:,:,:,end)),4)>=PDproc.WBMaskTh) = 1;
WMmask=zeros(size(squeeze(TPMs(:,:,:,1))));
WMmask(squeeze(TPMs(:,:,:,2))>=PDproc.WMMaskTh) = 1;

% Save masked A map for bias-field correction later
V_maskedA = spm_vol(fA);
V_maskedA.fname = fullfile(calcpath,['masked_' spm_str_manip(V_maskedA.fname,'t')]);
maskedA = spm_read_vols(spm_vol(fA)).*WBmask;
maskedA(maskedA==Inf) = 0;
maskedA(isnan(maskedA)) = 0;
maskedA(maskedA==threshA) = 0;
spm_write_vol(V_maskedA,maskedA);

% Bias-field correction of masked A map
clear matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel.vols = {V_maskedA.fname};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = PDproc.biasreg;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = PDproc.biasfwhm;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0];
output_list = spm_jobman('run', matlabbatch);

% Bias field correction of A map. 
% Bias field calculation is based on the masked A map, while correction
% must be applied to the unmasked A map. The BiasField is therefore
% retrieved from previous step and applied onto the original A map. 
BFfnam = output_list{1}.channel.biasfield{1};
BF = double(spm_read_vols(spm_vol(BFfnam)));
Y = BF.*spm_read_vols(spm_vol(fA));

% Calibration of flattened A map to % water content using typical white
% matter value from the litterature (69%)
A_WM = WMmask.*Y;
Y = Y/mean(A_WM(A_WM~=0))*69;
fprintf(1,'mean White Matter intensity: %.1f\n',mean(A_WM(A_WM~=0)));
fprintf(1,'SD White Matter intensity %.1f\n',std(A_WM(A_WM~=0),[],1));
Y(Y>200) = 0;
% MFC: Estimating Error for data set to catch bias field issues:
errorEstimate = std(A_WM(A_WM > 0))./mean(A_WM(A_WM > 0));
Vsave = spm_vol(fA);
Vsave.descrip = [Vsave.descrip '. Error Estimate: ', num2str(errorEstimate)];
if errorEstimate > 0.06 %#ok<BDSCI>
    % MFC: Testing on 15 subjects showed 6% is a good cut-off:
    warning(['Error estimate is high: ', Vsave.fname]);
end
if mpm_params.QA.enable
    if exist(mpm_params.QA.fnam,'file')
        mpm_params.QA = spm_jsonread(mpm_params.QA.fnam);
    end
    mpm_params.QA.PD.mean = mean(A_WM(A_WM > 0));
    mpm_params.QA.PD.SD = std(A_WM(A_WM > 0));
    spm_jsonwrite(mpm_params.QA.fnam,mpm_params.QA,struct('indent','\t'));
end
spm_write_vol(Vsave,Y);

end

%% =======================================================================%
% Sort out all parameters required for the MPM calculation. Check whether
% input data are coherent. Retrieves default parameters from the
% hmri_get_defaults. 
%=========================================================================%
function mpm_params = get_mpm_params(jobsubj)

% global parameters
mpm_params.json = hmri_get_defaults('json');
mpm_params.centre = hmri_get_defaults('centre');
mpm_params.calcpath = jobsubj.path.mpmpath;
mpm_params.respath = jobsubj.path.respath;
mpm_params.QA.enable = hmri_get_defaults('qMRI_maps.QA'); % quality assurance
if mpm_params.QA.enable
    mpm_params.QA.fnam = fullfile(mpm_params.respath,'MPM_map_creation_quality_assessment.json');
    spm_jsonwrite(mpm_params.QA.fnam,mpm_params.QA,struct('indent','\t'));
end
mpm_params.ACPCrealign = hmri_get_defaults('qMRI_maps.ACPCrealign'); % realigns qMRI maps to MNI

% retrieve input file names for map creation
mpm_params.input.MTw.fname   = char(jobsubj.raw_mpm.MT); % P_mtw
mpm_params.input.PDw.fname   = char(jobsubj.raw_mpm.PD); % P_pdw
mpm_params.input.T1w.fname   = char(jobsubj.raw_mpm.T1); % P_t1w

% maximum TE for averaging (ms)
mpm_params.input.TE_limit = 30; 

% acquisition parameters of MTw images
p = get_trtefa(mpm_params.input.MTw.fname);
if ~isempty(p)
    mpm_params.input.MTw.TE = cat(1,p.te);
    mpm_params.input.MTw.TR = p(1).tr;
    mpm_params.input.MTw.fa = p(1).fa;
else
    fprintf(1,'WARNING: No TE/TR/FA values found for MTw images. Fallback to defaults.\n');
    MPMacq = hmri_get_defaults('MPMacq');
    mpm_params.input.MTw.TE = MPMacq.TE_mtw;
    mpm_params.input.MTw.TR = MPMacq.TR_mtw;
    mpm_params.input.MTw.fa = MPMacq.fa_mtw;
end

% acquisition parameters of PDw images
p = get_trtefa(mpm_params.input.PDw.fname);
if ~isempty(p)
    mpm_params.input.PDw.TE = cat(1,p.te);
    mpm_params.input.PDw.TR = p(1).tr;
    mpm_params.input.PDw.fa = p(1).fa;
else
    fprintf(1,'WARNING: No TE/TR/FA values found for PDw images. Fallback to defaults.\n');
    MPMacq = hmri_get_defaults('MPMacq');
    mpm_params.input.MTw.TE = MPMacq.TE_pdw;
    mpm_params.input.MTw.TR = MPMacq.TR_pdw;
    mpm_params.input.MTw.fa = MPMacq.fa_pdw;
end

% acquisition parameters of T1w images
p = get_trtefa(mpm_params.input.T1w.fname);
if ~isempty(p)
    mpm_params.input.T1w.TE = cat(1,p.te);
    mpm_params.input.T1w.TR = p(1).tr;
    mpm_params.input.T1w.fa = p(1).fa;
else
    fprintf(1,'WARNING: No TE/TR/FA values found for T1w images. Fallback to defaults.\n');
    MPMacq = hmri_get_defaults('MPMacq');
    mpm_params.input.MTw.TE = MPMacq.TE_t1w;
    mpm_params.input.MTw.TR = MPMacq.TR_t1w;
    mpm_params.input.MTw.fa = MPMacq.fa_t1w;
end

% identify the protocol to define RF spoiling correction parameters
% retrieve all available protocols:
MPMacq_sets = hmri_get_defaults('MPMacq_set');
% current protocol is defined by [TR_pdw TR_t1w fa_pdw fa_t1w]:
MPMacq_prot = [mpm_params.input.PDw.TR mpm_params.input.T1w.TR mpm_params.input.PDw.fa mpm_params.input.T1w.fa]; 
% then match the values and find protocol tag
nsets = numel(MPMacq_sets.vals);
ii = 0; mtch = false;
while ~mtch && ii < nsets
    ii = ii+1;
    if all(MPMacq_prot == MPMacq_sets.vals{ii})
        mtch  = true;
        prot_tag = MPMacq_sets.tags{ii};
        fprintf(1,'INFO: MPM acquisition protocol = %s.\n', prot_tag);
    end
end
if ~mtch
    prot_tag = 'Unknown';
    fprintf(1,'WARNING: MPM protocol unknown. No RF spoiling correction will be applied.\n');
end
% now retrieve RF spoiling correction parameters
mpm_params.proc.RFC = hmri_get_defaults(['rfcorr.',prot_tag]);

% other processing parameters from defaults
% load threshold to save qMRI maps
mpm_params.proc.threshall = hmri_get_defaults('qMRI_maps_thresh');
% load PD maps processing parameters
mpm_params.proc.PD = hmri_get_defaults('PDproc');
% whether OLS R2* is calculated
mpm_params.proc.R2sOLS = hmri_get_defaults('R2sOLS');

% check that echo times are identical (common echoes only)
TE_mtw = mpm_params.input.MTw.TE;
TE_pdw = mpm_params.input.PDw.TE;
TE_t1w = mpm_params.input.T1w.TE;
nr_c_echoes = min([length(TE_mtw), length(TE_pdw), length(TE_t1w)]);
for nr = 1:nr_c_echoes
    if ~((TE_mtw(nr) == TE_pdw(nr)) && (TE_pdw(nr) == TE_t1w(nr)))
        error('Echo times do not match! Aborting.');
    end
end

end


%% =======================================================================%
% To arrange the metadata structure for MPM calculation output.
%=========================================================================%
function metastruc = init_mpm_output_metadata(input_files, mpm_params)

proc.descrip = 'MPM calculation';
proc.version = hmri_get_version;
proc.params = mpm_params;

% must be defined on the spot, default values here
output.imtype = 'MPM';
output.units = 'a.u.';

metastruc = init_output_metadata_structure(input_files, proc, output);

end


%% =======================================================================%
% To extract the TR/TE/FA values out of the descriptor field of 
% the nifti header or the extended header/JSON metadata.
%=========================================================================%
function p = get_trtefa(P)

N = nifti(P);
nN = numel(N);
p(nN) = struct('tr',[],'te',[],'fa',[]);

for ii = 1:numel(N)
    p(ii).tr = get_metadata_val(P(ii,:),'RepetitionTime');
    p(ii).te = get_metadata_val(P(ii,:),'EchoTime');
    p(ii).fa = get_metadata_val(P(ii,:),'FlipAngle');
end

end