function [fR1, fR2s, fMT, fA, PPDw, PT1w, PMTw]  = hmri_create_MTProt(jobsubj, P_trans, P_receiv) %#ok<*STOUT>
%==========================================================================
% This is hmri_create_MTProt, part of the hMRI-Toolbox.
%
% PURPOSE
% Estimation of quantitative parameters (R1, R2*, PD, MT) from
% multi-contrast multi-echo FLASH protocol 
% 
% FORMAT
% [fR1, fR2s, fMT, fA, PPDw, PT1w, PMTw]  = hmri_create_MTProt(jobsubj, P_trans, P_receiv)
%
% INPUTS
%   jobsubj     parameters for one subject out of the job list.
%               NB: ONE SINGLE DATA SET FROM ONE SINGLE SUBJECT IS
%               PROCESSED HERE, LOOP OVER SUBJECTS DONE AT HIGHER LEVEL.
%   P_trans     filenames of a pair of magnitude image + transmit bias map
%               [p.u.] of B1+ field (not mandatory) 
%   P_receive   filenames of a pair of magnitude image + sensitivity map of
%               phased-array coil relative to BC [p.u.] (not mandatory) 
%
% OUTPUTS
%   fR1     R1 map output filename (only quantitative if B1+ map provided)
%   fR2s    R2* map output filename (simple fit on PD data or OLS across
%           all contrasts, according to defaults settings)
%   fMT     Magnetization transfer (MT) map output filename  
%   fA      Proton density map output filename (free water concentration
%           (PD) [%] or signal amplitude (A) [a.u.] depending on defaults
%           settings (PDproc.PDmap)). 
%   PPDw    averate PD-weighted image filename (or OLS fit at TE=0 if fullOLS = true) 
%   PT1w    average T1-weighted image filename (or OLS fit at TE=0 if fullOLS = true)  
%   PMTw    average MT-weighted image filename (or OLS fit at TE=0 if fullOLS = true)  
%
% OTHER USEFUL VARIABLES EXPLAINED
%   P_mtw, P_pdw, P_t1w (from jobsubj.raw_mpm) are MTw, PDw, T1w series of
%           images respectively (multiple echoes) 
%   TE_mtw, TE_pdw, TE_t1w: echo times of the first echo (volume) of each
%           contrast, respectively.
%   TR_mtw, TR_pdw, TR_t1w: repetition time for each contrast,
%           respectively. 
%   fa_mtw, fa_pdw, fa_t1w: excitation flip angles for each contrast,
%           respectively. 
%
%==========================================================================
% FEATURES AND REFERENCES
%   Estimation of R1 and MT maps
%       Helms et al., Magnetic Resonance in Medicine 60:1396-1407 (2008)
%       Helms et al., Magnetic Resonance in Medicine 59:667-672 (2008)
%
%   B1 correction of MT maps
%       Weiskopf et al., Front. Neurosci. 2013, doi:
%       10.3389/fnins.2013.00095 
%
%   Imperfect RF spoiling correction
%       Preibisch and Deichmann, MRM 61:125-135 (2009)
%       Corrects for residual transverse magnetization coherences that
%       remain despite RF spoiling and result in deviations from the Ernst
%       equation (ideal expected signal in a FLASH acquisition). 
%       Modelling data and P2_a, P2_b correction factors calculation based
%       on code supplied by Ralf Deichmann. Only a small subset of
%       experimental parameters have been used and RF spoiling correction
%       can only be applied if these correction factors are defined
%       specifically for the acquisition protocol used.
%
%   OLS R2* map calculation 
%       Ordinary least square fit of R2* estimated from all echoes from all
%       three contrasts instead of the PD-weighted image only. This
%       increases the SNR in R2* maps. It is noted that it potentially
%       mixes in additional MT and T1 contrast, as discussed in: 
%       Weiskopf et al. Front. Neurosci. 2014 DOI: 10.3389/fnins.2014.00278 
%       This reference should be cited if you use this output.
%
%==========================================================================
% Written by
%   Gunther Helms, MR-Research in Neurology and Psychiatry, University
%       of Goettingen 
%   Martina Callaghan, John Ashburner, Wellcome Trust Centre for
%       Neuroimaging at UCL, London  
%   Antoine Lutti, CHUV, Lausanne, Switzerland
%   Nikolaus Weiskopf, Department of Neurophysics, Max Planck Institute 
%       for Human Cognitive and Brain Sciences, Leipzig, Germany 
%
%==========================================================================

%% =======================================================================%
% Initiate processing
%=========================================================================%
fprintf(1, ['---------------------- MPM PROCESSING ----------------------\n' ...
    'Create maps from multi-contrast multi-echo FLASH protocol.\n']); 
                    
% retrieves all required parameters for MPM processing
mpm_params = get_mpm_params(jobsubj);

% for convenience, define a few parameters to make formulae more readable
% and avoid number of repetitions:

% index number for each contrast - zero index means no images available
PDidx = mpm_params.PDidx;
T1idx = mpm_params.T1idx;
MTidx = mpm_params.MTidx;
% TE/TR/FA for each contrast 
% T1w and MTw contrasts are optional.
% PDw must be present or nothing can be calculated. The script will abort
% at the next line if no PDw echoes available. In that case, a warning has
% been thrown earlier (in get_mpm_params).
TE_pdw = mpm_params.input(PDidx).TE;
TR_pdw = mpm_params.input(PDidx).TR;
fa_pdw = mpm_params.input(PDidx).fa;
if T1idx
    TE_t1w = mpm_params.input(T1idx).TE; %#ok<NASGU>
    TR_t1w = mpm_params.input(T1idx).TR;
    fa_t1w = mpm_params.input(T1idx).fa;
end
if MTidx
    TE_mtw = mpm_params.input(MTidx).TE; %#ok<NASGU>
    TR_mtw = mpm_params.input(MTidx).TR;
    fa_mtw = mpm_params.input(MTidx).fa;
end
% other parameters
threshall = mpm_params.proc.threshall;
PDproc = mpm_params.proc.PD;
RFC = mpm_params.proc.RFC;
dt = [spm_type('float32'),spm_platform('bigend')]; % for nifti output
outbasename = spm_file(mpm_params.input(end).fnam(1,:),'basename'); % for all output files
calcpath = mpm_params.calcpath;
mpm_params.outbasename = outbasename;
respath = mpm_params.respath;
supplpath = mpm_params.supplpath;
% Number of echoes averaged (maximum number or echoes available for ALL
% contrasts AND under TE_limit (+1) - see get_mpm_params)
avg_nr = mpm_params.nr_echoes4avg; 

%% =======================================================================%
% Calculate R2* map from PDw echoes
%=========================================================================%
fprintf(1,'\n    -------- R2* map calculation --------\n');

% load PDw images
V_pdw = spm_vol(mpm_params.input(PDidx).fnam);
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
        Y  = Y + W(i)*log(max(spm_slice_vol(V_pdw(i),M1,dm(1:2),mpm_params.interp),1));
    end
    Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

% Set and write metadata
input_files = mpm_params.input(PDidx).fnam;
Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
Output_hdr.history.output.imtype = 'R2* map';
Output_hdr.history.output.units = 'ms-1';
set_metadata(fR2s,Output_hdr,mpm_params.json);


%% =======================================================================%
% Reading and averaging the images 
%=========================================================================%
fprintf(1,'\n    -------- Reading and averaging the images --------\n');

% Average only first few echoes for increased SNR and fit T2* 
Pavg = cell(1,mpm_params.ncon);
for ccon=1:mpm_params.ncon % loop over available contrasts
    Pavg{ccon}  = fullfile(calcpath,[outbasename '_' mpm_params.input(ccon).tag 'w.nii']);
    V           = spm_vol(mpm_params.input(ccon).fnam);
    dm          = V(1).dim;
    Ni          = nifti;
    Ni.mat      = V(1).mat;
    Ni.mat0     = V(1).mat;
    Ni.descrip  = sprintf('Averaged %sw images', mpm_params.input(ccon).tag);
    Ni.dat      = file_array(Pavg{ccon},dm,dt,0,1,0);
    create(Ni);
    spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
    % sm = 0;
    for p = 1:dm(3)
        M = spm_matrix([0 0 p]);
        Y = zeros(dm(1:2));
        for nr = 1:avg_nr
            M1 = V(nr).mat\V(1).mat*M;
            Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),mpm_params.interp);
        end
        Ni.dat(:,:,p) = Y/avg_nr;
        % sm = sm + sum(Y(:))/avg_nr;
        spm_progress_bar('Set',p);
    end
    % avg = sm/prod(dm);
    spm_progress_bar('Clear');
    
    input_files = mpm_params.input(ccon).fnam;
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = Ni.descrip;
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(Pavg{ccon},Output_hdr,mpm_params.json);
end

% Average T1w image for PD calculation 
% (average over PDproc.nr_echoes_forA echoes, see hmri_defaults):
if (PDidx && T1idx)
    PT1w_forA = fullfile(calcpath,[outbasename '_T1w_forA.nii']);
    V           = spm_vol(mpm_params.input(T1idx).fnam);
    dm          = V(1).dim;
    Ni          = nifti;
    Ni.mat      = V(1).mat;
    Ni.mat0     = V(1).mat;
    Ni.descrip  = 'Averaged T1w images for PD calculation';
    Ni.dat      = file_array(PT1w_forA,dm,dt, 0,1,0);
    create(Ni);
    spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
    for p = 1:dm(3),
        M = spm_matrix([0 0 p]);
        Y = zeros(dm(1:2));
        for nr = 1:PDproc.nr_echoes_forA,
            M1 = V(nr).mat\V(1).mat*M;
            Y  = Y + spm_slice_vol(V(nr),M1,dm(1:2),mpm_params.interp);
        end
        Ni.dat(:,:,p) = Y/PDproc.nr_echoes_forA;
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');
end

%% =======================================================================%
% Coregistering the images
%=========================================================================%
fprintf(1,'\n    -------- Coregistering the images  --------\n');

x_MT2PD = [];
if MTidx; x_MT2PD = coreg_mt(Pavg{PDidx}, Pavg{MTidx}); end
x_T12PD = [];   
if T1idx; 
    x_T12PD = coreg_mt(Pavg{PDidx}, Pavg{T1idx});
    coreg_mt(Pavg{PDidx}, PT1w_forA);
end

V_trans = [];
if ~isempty(P_trans)
    % Load B1 mapping data if available and coregister to PDw
    % P_trans(1,:) = magnitude image (anatomical reference for coregistration) 
    % P_trans(2,:) = B1 map (p.u.)
    coreg_bias_map(Pavg{PDidx}, P_trans);
    V_trans = spm_vol(P_trans);
end

V_receiv   = [];
if ~isempty(P_receiv)
    % Load sensitivity map if available and coregister to PDw
    % P_receiv(1,:) = magnitude image (anatomical reference for coregistration) 
    % P_receiv(2,:) = sensitivity map
    coreg_bias_map(Pavg{PDidx}, P_receiv);
    V_receiv = spm_vol(P_receiv);
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
Vavg(PDidx) = spm_vol(Pavg{PDidx});
if MTidx; Vavg(MTidx) = spm_vol(Pavg{MTidx}); end
if T1idx
    Vavg(T1idx) = spm_vol(Pavg{T1idx}); 
    VT1w_forA = spm_vol(PT1w_forA);
end


%% =======================================================================%
% multi-contrast R2* map calculation for quality assessment by AL
% METHOD: R2s map calculated for each contrast separately (PDw, MTw, T1w)
% to evaluate the SD of each R2* map in white matter (i.e. intra-run motion
% measurement).
%=========================================================================%
if mpm_params.QA.enable
    fprintf(1,'\n    -------- multi-contrast R2* map calculation for QA --------\n');
    
    for ccon = 1:mpm_params.ncon
        dt        = [spm_type('float32'),spm_platform('bigend')];
        Ni        = nifti;
        Ni.mat    = V_pdw(1).mat;
        Ni.mat0   = V_pdw(1).mat;
        Ni.descrip='OLS R2* map [1/ms]';
        Ni.dat    = file_array(fullfile(calcpath,[outbasename '_R2s_' mpm_params.input(ccon).tag 'w.nii']),dm,dt, 0,1,0);
        create(Ni);
        
        TE = mpm_params.input(ccon).TE;
        Vcon = spm_vol(mpm_params.input(ccon).fnam);
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
            
            for cecho = 1:size(TE,1)
                % Take slice p (defined in M) and map to a location in the
                % appropriate contrast using the matField entry for that
                % contrast, which has been co-registered to the PD-weighted
                % data:
                M1 = Vavg(ccon).mat\V_pdw(1).mat*M;
                
                % Third order B-spline interpolation for OLS R2* estimation
                % since we no longer assume that the echoes are perfectly
                % aligned as we do for the standard PDw derived R2* estimate.
                data(cecho,:,:) = log(max(spm_slice_vol(Vcon(cecho),M1,dm(1:2),mpm_params.interp),1));
            end
            Y = W*reshape(data, [size(TE,1) prod(dm(1:2))]);
            Y = -reshape(Y(end,:), dm(1:2));
            
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
        
    % OLS fit at TE=0: to be used instead of averaged echoes of each
    % contrast if "fullOLS" option is enabled
    if mpm_params.fullOLS
        fprintf(1,'\n    -------- and fit to TE=0 for all contrasts --------\n');

        Nmap = nifti;
        Pte0 = cell(1,mpm_params.ncon);
        for ccon = 1:mpm_params.ncon
            Pte0{ccon}  = fullfile(calcpath,[outbasename '_' mpm_params.input(ccon).tag 'w_OLSfit_TEzero.nii']);
            Ni          = nifti;
            Ni.mat      = V_pdw(1).mat;
            Ni.mat0     = V_pdw(1).mat;
            Ni.descrip  = sprintf('OLS fit to TE=0 for %sw images', mpm_params.input(ccon).tag);
            Ni.dat      = file_array(Pte0{ccon},dm,dt,0,1,0);
            create(Ni);
            
            % set metadata
            input_files = mpm_params.input(ccon).fnam;
            Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
            Output_hdr.history.output.imtype = Ni.descrip;
            Output_hdr.history.output.units = 'a.u.';
            set_metadata(Pte0{ccon},Output_hdr,mpm_params.json);            
            
            % re-load the updated NIFTI file (in case extended header has
            % been added, the offset has changed and must be updated before
            % writing the data to the file!)
            Nmap(ccon) = nifti(Pte0{ccon});
        end
    end % init nifti objects for fullOLS case
    
    fR2s_OLS    = fullfile(calcpath,[outbasename '_R2s_OLS' '.nii']);
    Ni          = nifti;
    Ni.mat      = V_pdw(1).mat;
    Ni.mat0     = V_pdw(1).mat;
    Ni.descrip  = 'OLS R2* map [1/ms]';
    Ni.dat      = file_array(fR2s_OLS,dm,dt,0,1,0);
    create(Ni);
    
    % Combine the data and echo times:
    nechoes = zeros(1,mpm_params.ncon);
    for ccon = 1:mpm_params.ncon
        nechoes(ccon) = size(mpm_params.input(ccon).fnam,1);
    end
    
    % Same formalism as for PDw fit but now extra colums for the "S(0)" 
    % amplitudes of the different contrasts:
    reg = zeros(sum(nechoes),mpm_params.ncon+1);
    for ccon = 1:mpm_params.ncon
        reg(sum(nechoes(1:ccon-1))+(1:nechoes(ccon)),ccon) = 1;
        reg(sum(nechoes(1:ccon-1))+(1:nechoes(ccon)),end) = mpm_params.input(ccon).TE;
    end
    W = (reg'*reg)\reg';
    
    spm_progress_bar('Init',dm(3),'OLS R2* fit','planes completed');
    for p = 1:dm(3),
        M = spm_matrix([0 0 p 0 0 0 1 1 1]);
        data = zeros([sum(nechoes) dm(1:2)]);
        
        for ccon = 1:mpm_params.ncon
            
            Vcon = spm_vol(mpm_params.input(ccon).fnam);
            
            % Take slice p (defined in M) and map to a location in the
            % appropriate contrast using the V.mat field entry for that
            % contrast, which has been co-registered to the PD-weighted
            % data:
            M1 = Vavg(ccon).mat\V_pdw(1).mat*M;
            
            for cecho = 1:nechoes(ccon)               
                % Third order B-spline interpolation for OLS R2* estimation
                % since we no longer assume that the echoes are perfectly
                % aligned as we do for the standard PDw derived R2*
                % estimate. 
                data(sum(nechoes(1:ccon-1))+cecho,:,:) = log(max(spm_slice_vol(Vcon(cecho),M1,dm(1:2),mpm_params.interp),1));
            end
        end
        Y = W*reshape(data, [sum(nechoes) prod(dm(1:2))]);

        % Writes "fullOLS" images (OLS fit to TE=0 for each contrast)
        if mpm_params.fullOLS
            for ccon = 1:mpm_params.ncon
                Nmap(ccon).dat(:,:,p) = reshape(exp(Y(ccon,:)), dm(1:2));
            end
        end
        
        Y = -reshape(Y(end,:), dm(1:2));
        
        % NB: mat field defined by V_pdw => first PDw echo
        % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), 
        % negative values are allowed to preserve Gaussian distribution.
        Ni.dat(:,:,p) = max(min(Y,threshall.R2s),-threshall.R2s); 
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');

    % Set metadata (R2S_OLS)
    input_files = mpm_params.input(PDidx).fnam;
    if (T1idx); input_files = char(input_files, mpm_params.input(T1idx).fnam); end
    if (MTidx); input_files = char(input_files, mpm_params.input(MTidx).fnam); end
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = 'R2*-OLS map';
    Output_hdr.history.output.units = 'ms-1';
    set_metadata(fullfile(calcpath,[outbasename '_R2s_OLS' '.nii']),Output_hdr,mpm_params.json);
   
end % OLS code

% if "fullOLS" option enabled, Pte0 images replace Pavg images in the rest
% of the code. We need to apply the substitution and reload the images...
if mpm_params.fullOLS
    for ccon = 1:mpm_params.ncon
        Pavg{ccon} = Pte0{ccon};
        Vavg(ccon) = spm_vol(Pavg{ccon});
    end
end
        

%% =======================================================================%
% Prepare output for R1, PD and MT maps
%=========================================================================%
% description fields and file names of output images

coutput = 0;

if T1idx
    coutput = coutput+1;
    output_suffix{coutput} = 'R1';
    units{coutput} = '1000/s';
    if isempty(V_trans)
        descrip{coutput} = 'R1 map (no B1+ bias correction applied)';
    else
        descrip{coutput} = 'R1 map (with B1+ bias correction)';
    end
    R1map_idx = coutput;
end

if T1idx
    coutput = coutput+1;
    if PDproc.PDmap && ~isempty(V_trans)
        output_suffix{coutput} = 'PD';
        descrip{coutput} = 'Water concentration [%]';
        units{coutput} = '%';
    else
        output_suffix{coutput} = 'A';
        descrip{coutput} = 'Signal amplitude [a.u.]';
        units{coutput} = 'a.u.';
    end
    Amap_idx = coutput;
end

if (MTidx && T1idx)
    coutput = coutput+1;
    output_suffix{coutput} = 'MT';
    descrip{coutput} = 'Delta MT map';
    units{coutput} = 'a.u.';
    MTmap_idx = coutput;
end

if (MTidx && PDidx)
    if (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw) % additional MTR image...
        coutput = coutput+1;
        output_suffix{coutput} = 'MTR';
        descrip{coutput} = 'Classic MTR image';
        units{coutput} = 'a.u.';
        MTRmap_idx = coutput;    
    end
end

noutput = coutput;

% define NIFTI objects for output images
Nmap    = nifti;
for ii=1:noutput
    dm         = V_pdw(1).dim;
    Ni         = nifti;
    Ni.mat     = V_pdw(1).mat;
    Ni.mat0    = V_pdw(1).mat;
    Ni.descrip = descrip{ii};
    Ni.dat     = file_array(fullfile(calcpath,[outbasename '_' output_suffix{ii} '.nii']),dm,dt, 0,1,0);
    create(Ni);
    Nmap(ii) = Ni;
end

fR1 = '';
fA = '';
fMT = '';
if (PDidx && T1idx)
    fR1 = fullfile(calcpath,[outbasename '_' output_suffix{R1map_idx} '.nii']);
    fA  = fullfile(calcpath,[outbasename '_' output_suffix{Amap_idx} '.nii']);
    if MTidx
        fMT = fullfile(calcpath,[outbasename '_' output_suffix{MTmap_idx} '.nii']);
    end        
end


%% =======================================================================%
% Map calculation continued (R1, PD, MT) 
%=========================================================================%
fprintf(1,'\n    -------- Map calculation continued (R1, PD, MT) --------\n');

M0 = Ni.mat;
dm = size(Ni.dat);

fa_pdw_rad = fa_pdw * pi / 180;
if MTidx; fa_mtw_rad = fa_mtw * pi / 180; end
if T1idx; fa_t1w_rad = fa_t1w * pi / 180; end

spm_progress_bar('Init',dm(3),'Calculating maps','planes completed');

for p = 1:dm(3)
    M = M0*spm_matrix([0 0 p]);

    % PDw images are always available, so this bit is always loaded:
    PDw = spm_slice_vol(Vavg(PDidx),Vavg(PDidx).mat\M,dm(1:2),mpm_params.interp);
    
    if ~isempty(V_trans)
        f_T = spm_slice_vol(V_trans(2,:),V_trans(2,:).mat\M,dm(1:2),mpm_params.interp)/100; % divide by 100, since p.u. maps
    else
        f_T = [];
    end
    if ~isempty(V_receiv) && ~isempty(V_trans)
        f_R = spm_slice_vol(V_receiv(2,:),V_receiv(2,:).mat\M,dm(1:2),mpm_params.interp)/100; % divide by 100, since p.u. maps
        f_R = f_R .* f_T; % f_R is only the sensitivity map and not the true receive bias map, therefore needs to be multiplied by transmit bias (B1+ approx. B1- map)
    else
        f_R = [];
    end
    
    % Standard magnetization transfer ratio (MTR) in percent units [p.u.]
    % only if trpd = trmt and fapd = famt and if PDw and MTw images are
    % available
    if (MTidx && PDidx)
        MTw = spm_slice_vol(Vavg(MTidx),Vavg(MTidx).mat\M,dm(1:2),mpm_params.interp);
        if (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw) % additional MTR image...
            MTR = (PDw-MTw)./(PDw+eps) * 100;
            % write MTR image
            Nmap(MTRmap_idx).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
        end          
    end
    
    % T1 map and A/PD maps can only be calculated if T1w images are
    % available:
    if T1idx

        T1w = spm_slice_vol(Vavg(T1idx),Vavg(T1idx).mat\M,dm(1:2),mpm_params.interp);
        
        % if "fullOLS" option enabled, use the OLS fit at TE=0 as
        % "T1w_forA"; otherwise use the average calculated earlier (by
        % default, corresponds to the first echo to reduce R2* bias)
        if mpm_params.fullOLS
            T1w_forA = T1w;
        else
            T1w_forA = spm_slice_vol(VT1w_forA,VT1w_forA.mat\M,dm(1:2),mpm_params.interp);
        end
    
        if isempty(f_T)
            % semi-quantitative T1
            T1 = ((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)) ./ ...
                max((T1w * (fa_t1w_rad / 2 / TR_t1w)) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps);
            R1 = (((T1w * (fa_t1w_rad / 2 / TR_t1w)) - (PDw * fa_pdw_rad / 2 / TR_pdw)) ./ ...
                max(((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)),eps))*10^6;
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
                    ((((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)+eps) ./ ...
                    max((T1w * fa_t1w_rad / 2 / TR_t1w) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps))./f_T.^2);
            else
                % MFC: We do not have P2_a or P2_b parameters for this sequence
                % => T1 = T1app
                T1 = ((((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)+eps) ./ ...
                    max((T1w * fa_t1w_rad / 2 / TR_t1w) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps))./f_T.^2);
            end
            
            R1 = 1./T1*10^6;
        end
        
        T1       = max(T1,0);
        R1(R1<0) = 0;
        tmp      = R1;
        Nmap(R1map_idx).dat(:,:,p) = min(max(tmp,-threshall.R1),threshall.R1); % truncating images
        
        % A values proportional to PD
        if (~isempty(f_T)) && (~isempty(f_R))
            % Transmit and receive bias corrected quantitative A values
            % again: correct A for transmit bias f_T and receive bias f_R
            % Acorr = A / f_T / f_R , proportional PD
            A = (T1 .* (T1w_forA * fa_t1w_rad / 2 / TR_t1w) + (T1w_forA / fa_t1w_rad))./f_T./f_R;
        elseif(~isempty(f_T))&&(isempty(f_R))%&&(PDproc.PDmap)
            A = T1 .* (T1w_forA .*(fa_t1w_rad*f_T) / 2 / TR_t1w) + (T1w_forA ./ (fa_t1w_rad*f_T));
        else
            % semi-quantitative A
            A = T1 .* (T1w_forA * fa_t1w_rad / 2 / TR_t1w) + (T1w_forA / fa_t1w_rad);
        end
        
        tmp      = A;
        Nmap(Amap_idx).dat(:,:,p) = max(min(tmp,threshall.A),-threshall.A);
        % dynamic range increased to 10^5 to accommodate phased-array coils and symmetrical for noise distribution
        
        % for MT maps calculation, one needs MTw images on top of the T1w
        % and PDw ones...
        if MTidx
            MTw = spm_slice_vol(Vavg(MTidx),Vavg(MTidx).mat\M,dm(1:2),3);
            T1_forMT = ((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)) ./ ...
                max((T1w * (fa_t1w_rad / 2 / TR_t1w)) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps);
            A_forMT = T1_forMT .* (T1w * fa_t1w_rad / 2 / TR_t1w) + (T1w / fa_t1w_rad);
            
            % MT in [p.u.]; offset by - famt * famt / 2 * 100 where MT_w = 0 (outside mask)
            MT       = ( (A_forMT * fa_mtw_rad - MTw) ./ (MTw+eps) ./ (T1_forMT + eps) * TR_mtw - fa_mtw_rad^2 / 2 ) * 100;
            if (~isempty(f_T))
                MT = MT .* (1 - 0.4) ./ (1 - 0.4 * f_T);
            end
            
            tmp      = MT;
            Nmap(MTmap_idx).dat(:,:,p) = max(min(tmp,threshall.MT),-threshall.MT);
        end
    end
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

% set metadata for all output images
input_files = mpm_params.input(PDidx).fnam;
if (T1idx); input_files = char(input_files, mpm_params.input(T1idx).fnam); end
if (MTidx); input_files = char(input_files, mpm_params.input(MTidx).fnam); end
Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
for ctr = 1:noutput
    Output_hdr.history.output.imtype = descrip(ctr);
    Output_hdr.history.output.units = units(ctr);
    set_metadata(fullfile(calcpath,[outbasename '_' output_suffix{ctr} '.nii']),Output_hdr,mpm_params.json);
end

%% =======================================================================%
% ACPC Realign all images - only if MT map created
%=========================================================================%
if mpm_params.ACPCrealign 
    if (MTidx && PDidx && T1idx)
    fprintf(1,'\n    -------- ACPC Realign all images --------\n');
    
    % Define and calculate masked MT image
    % Load MT image
    V_MT = spm_vol(fMT);
    MTimage = spm_read_vols(V_MT);    
    % Define new file name for masked MT image
    V_MT.fname = fullfile(calcpath,['masked_' spm_str_manip(fMT,'t')]);
    % Load average PDw image (mask based on averaged PDw image)
    PDWimage = spm_read_vols(Vavg(PDidx));
    % Mask MT image and save the masked MT image ('masked_..._MT.nii')
    MTimage(PDWimage<0.6*mean(PDWimage(:)))=0;
    spm_write_vol(V_MT,MTimage);

    % Use masked MT image to calculate transformation for ACPC realignment
    % (to increase robustness in segmentation):
    [~,R] = hmri_create_comm_adjust(1,V_MT.fname,V_MT.fname,8,0,fullfile(spm('Dir'),'canonical','avg152T1.nii')); 
    
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
    spm_jsonwrite(fullfile(supplpath,'MPM_map_creation_ACPCrealign_transformation_matrix.json'),R,struct('indent','\t'));

    else
        fprintf(1,['\nWARNING: ACPC Realign was enabled, but no MT map was available \n' ...
                   'to proceed. ACPC realignment must be done separately, e.g. you can \n'...
                   'run [hMRI tools > Auto-reorient] before calculating the maps.\n' ...
                   'NOTE: segmentation might crash if no initial reorientation.\n']);
    end
end

% for quality assessment and/or PD map calculation
% segmentation preferentially performed on MT map but can be done on R1 map
% if no MT map available. Therefore, we must at least have R1 available,
% i.e. both PDw and T1w inputs...
if (mpm_params.QA.enable||(PDproc.PDmap)) && (PDidx && T1idx)
    if ~isempty(fMT); 
        Vsave = spm_vol(fMT);
    else % ~isempty(fR1); 
        Vsave = spm_vol(fR1); 
    end
    MTtemp = spm_read_vols(Vsave);
    %The 5 outer voxels in all directions are nulled in order to remove artefactual effects from the MT map on segmentation:
    MTtemp(1:5,:,:)=0; MTtemp(end-5:end,:,:)=0;
    MTtemp(:,1:5,:)=0; MTtemp(:,end-5:end,:)=0;
    MTtemp(:,:,1:5)=0; MTtemp(:,:,end-5:end)=0;
    Vsave.fname = spm_file(Vsave.fname,'suffix','_outer_suppressed');
    spm_write_vol(Vsave,MTtemp);
    
    % all segmentation steps should now use US and identical TPMs for
    % coherence and uniformity across the toolbox...
    clear matlabbatch
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {Vsave.fname};
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    output_list = spm_jobman('run', matlabbatch);
    fTPM = char(cat(1,output_list{1}.tiss.c));
end

% for quality assessment - the above segmentation must have run
if mpm_params.QA.enable && exist('fTPM','var')
    
    % Load existing QA results
    if exist(mpm_params.QA.fnam,'file')
        mpm_params.QA = spm_jsonread(mpm_params.QA.fnam);
    end

    % Calculate WM mask
    TPMs = spm_read_vols(spm_vol(fTPM));
    WMmask = zeros(size(squeeze(TPMs(:,:,:,2))));
    WMmask(squeeze(TPMs(:,:,:,2))>=PDproc.WMMaskTh) = 1;        
    WMmask = spm_erode(spm_erode(double(WMmask)));
    
    % Load OLS R2s maps calculated for each contrast, mask them and
    % calculate SD within the WM mask (measure of the intra-run motion for
    % each contrast)
    for ccon = 1:mpm_params.ncon
        R2s = spm_read_vols(spm_vol(spm_select('FPList',calcpath,sprintf('^s.*_R2s_%sw.nii$',mpm_params.input(ccon).tag))));
        MaskedR2s = squeeze(R2s.*WMmask);
        SDR2s = std(MaskedR2s(MaskedR2s~=0),[],1);
        mpm_params.QA.SDR2s.([mpm_params.input(ccon).tag 'w']) = SDR2s;
    end
    
    spm_jsonwrite(mpm_params.QA.fnam, mpm_params.QA, struct('indent','\t'));
end

% PD map calculation
if ~isempty(f_T) && isempty(f_R) && PDproc.PDmap
    PDcalculation(fA,fTPM,mpm_params);
end

% copy final result files into Results directory
% NB: to avoid ambiguity for users, only the 4 final maps to be used in
% further analysis are in Results, all other files (additional maps, 
% processing parameters, etc) are in Results/Supplementary.
if ~isempty(fR1)
    fR1_final = fullfile(respath, spm_file(fR1,'filename'));
    copyfile(fR1,fR1_final);
    try copyfile([spm_str_manip(fR1,'r') '.json'],[spm_str_manip(fR1_final,'r') '.json']); end %#ok<*TRYNC>
    fR1 = fR1_final;
end

fR2s_final = fullfile(respath, spm_file(fR2s,'filename'));
copyfile(fR2s,fR2s_final);
try copyfile([spm_str_manip(fR2s,'r') '.json'],[spm_str_manip(fR2s_final,'r') '.json']); end
fR2s = fR2s_final;

% NB: if OLS calculation of R2s map has been done, the output file for R2s
% map is the OLS result. In that case, the basic R2s map is moved to
% Results/Supplementary while the R2s_OLS is copied into results directory: 
if mpm_params.proc.R2sOLS
    % move basin R2s map to Results/Supplementary
    movefile(fR2s_final, fullfile(supplpath, spm_file(fR2s_final,'filename')));
    try movefile([spm_str_manip(fR2s_final,'r') '.json'],fullfile(supplpath, [spm_file(fR2s_final,'basename') '.json'])); end
    % copy OLS_R2s map to Results
    fR2s_OLS_final = fullfile(respath, spm_file(fR2s_OLS,'filename'));
    copyfile(fR2s_OLS,fR2s_OLS_final);
    try copyfile([spm_str_manip(fR2s_OLS,'r') '.json'],[spm_str_manip(fR2s_OLS_final,'r') '.json']); end
    % the hmri_create_MTProt fR2s output in now the R2s_OLS map
    fR2s = fR2s_OLS_final;
end

if ~isempty(fMT)
    fMT_final = fullfile(respath, spm_file(fMT,'filename'));
    copyfile(fMT,fMT_final);
    try copyfile([spm_str_manip(fMT,'r') '.json'],[spm_str_manip(fMT_final,'r') '.json']); end
    fMT = fMT_final;
end

if ~isempty(fA)
    fA_final = fullfile(respath, spm_file(fA,'filename'));
    copyfile(fA,fA_final);
    try copyfile([spm_str_manip(fA,'r') '.json'],[spm_str_manip(fA_final,'r') '.json']); end
    fA = fA_final;
end

PPDw_final = fullfile(supplpath, spm_file(Pavg{PDidx},'filename'));
copyfile(Pavg{PDidx},PPDw_final);
try copyfile([spm_str_manip(Pavg{PDidx},'r') '.json'],[spm_str_manip(PPDw_final,'r') '.json']); end
PPDw = PPDw_final;

if T1idx
    PT1w_final = fullfile(supplpath, spm_file(Pavg{T1idx},'filename'));
    copyfile(Pavg{T1idx},PT1w_final);
    try copyfile([spm_str_manip(Pavg{T1idx},'r') '.json'],[spm_str_manip(PT1w_final,'r') '.json']); end
    PT1w = PT1w_final;
else
    PT1w = '';
end

if MTidx
    PMTw_final = fullfile(supplpath, spm_file(Pavg{MTidx},'filename'));
    copyfile(Pavg{MTidx},PMTw_final);
    try copyfile([spm_str_manip(Pavg{MTidx},'r') '.json'],[spm_str_manip(PMTw_final,'r') '.json']); end
    PMTw = PMTw_final;
else
    PMTw = '';
end

% save processing params (mpm_params)
spm_jsonwrite(fullfile(supplpath,'MPM_map_creation_mpm_params.json'),mpm_params,struct('indent','\t'));

spm_progress_bar('Clear');

end

%% =======================================================================%
% To coregister the structural images for MT protocol 
%=========================================================================%
function [x] = coreg_mt(P_ref, P_src)

for src_nr = 1:size(P_src,1)
    VG = spm_vol(P_ref);
    VF = spm_vol(P_src(src_nr,:));
    coregflags.sep = [4 2];
    x = spm_coreg(VG,VF, coregflags);
    M  = inv(spm_matrix(x));
    MM = spm_get_space(deblank(VF.fname));
    spm_get_space(deblank(deblank(VF.fname)), M*MM); %#ok<*MINV>
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
function PDcalculation(fA, fTPM, mpm_params)
% fA is the filename of the output A image 
% (not yet properly quantitative PD map when it enters PDcalculation)
% fTMP is the list of TPMs generated from MT map
fprintf(1,'\n    -------- Proton density map calculation --------\n');

PDproc = mpm_params.proc.PD;
threshA = mpm_params.proc.threshall.A;
calcpath = mpm_params.calcpath;

TPMs = spm_read_vols(spm_vol(fTPM));
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
fprintf(1,'\nINFO (PD calculation):\n\tmean White Matter intensity: %.1f\n',mean(A_WM(A_WM~=0)));
fprintf(1,'\tSD White Matter intensity %.1f\n\n',std(A_WM(A_WM~=0),[],1));
Y(Y>200) = 0;
% MFC: Estimating Error for data set to catch bias field issues:
errorEstimate = std(A_WM(A_WM > 0))./mean(A_WM(A_WM > 0));
Vsave = spm_vol(fA);
Vsave.descrip = [Vsave.descrip '. Error Estimate: ', num2str(errorEstimate)];
if errorEstimate > 0.06 %#ok<BDSCI>
    % MFC: Testing on 15 subjects showed 6% is a good cut-off:
    fprintf(1,['\nWARNING: Error estimate is high for calculated PD map:\n%s' ...
        '\nError higher than 6% may indicate motion.\n'], Vsave.fname);
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
mpm_params.supplpath = jobsubj.path.supplpath;
mpm_params.QA.enable = hmri_get_defaults('qMRI_maps.QA'); % quality assurance
if mpm_params.QA.enable
    mpm_params.QA.fnam = fullfile(mpm_params.supplpath,'MPM_map_creation_quality_assessment.json');
    spm_jsonwrite(mpm_params.QA.fnam,mpm_params.QA,struct('indent','\t'));
end
mpm_params.ACPCrealign = hmri_get_defaults('qMRI_maps.ACPCrealign'); % realigns qMRI maps to MNI
mpm_params.interp = hmri_get_defaults('interp');
mpm_params.fullOLS = hmri_get_defaults('fullOLS'); % uses all echoes for OLS fit at TE=0

% retrieve input file names for map creation.
% the "mpm_params.input" field is an array, each element corresponds to a
% contrast.  
% if no input files for a given contrast, no input entry created and
% warning is thrown.
ccon = 0;
fprintf(1,'\nINFO: FLASH echoes loaded for each contrast are: ');
% 1) try MTw contrast:
tmpfnam   = char(jobsubj.raw_mpm.MT); % P_mtw
if isempty(tmpfnam)
    fprintf(1,'\n\t- WARNING: no MT-weighted FLASH echoes available!');
    mpm_params.MTidx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    fprintf(1,'\n\t- MT-weighted: %d echoes', size(tmpfnam,1));
    mpm_params.input(ccon).fnam = tmpfnam;
    mpm_params.input(ccon).tag = 'MT';  
    mpm_params.MTidx = ccon;
end  
% 2) try PDw contrast:
tmpfnam   = char(jobsubj.raw_mpm.PD); % P_pdw
if isempty(tmpfnam)
    fprintf(1,'\n\n\t- WARNING: no PD-weighted FLASH echoes available! \n\t\tThe map creation won''t be able to proceed!\n');
    mpm_params.PDidx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    fprintf(1,'\n\t- PD-weighted: %d echoes', size(tmpfnam,1));
    mpm_params.input(ccon).fnam = tmpfnam;
    mpm_params.input(ccon).tag = 'PD';  
    mpm_params.PDidx = ccon;
end  
% 3) try T1w contrast:
tmpfnam   = char(jobsubj.raw_mpm.T1); % P_t1w
if isempty(tmpfnam)
    fprintf(1,'\n\t- WARNING: no T1-weighted FLASH echoes available!');
    mpm_params.T1idx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    fprintf(1,'\n\t- T1-weighted: %d echoes', size(tmpfnam,1));
    mpm_params.input(ccon).fnam = tmpfnam;
    mpm_params.input(ccon).tag = 'T1';  
    mpm_params.T1idx = ccon; % zero index means no contrast available    
end 
mpm_params.ncon = ccon; % number of contrasts available
fprintf(1,'\n');


% collect TE, TR and FA for each available contrast
for ccon = 1:mpm_params.ncon
    p = get_trtefa(mpm_params.input(ccon).fnam);
    if ~isempty(p)
        mpm_params.input(ccon).TE = cat(1,p.te);
        mpm_params.input(ccon).TR = p(1).tr;
        mpm_params.input(ccon).fa = p(1).fa;
    else
        fprintf(1,'\nWARNING: No TE/TR/FA values found for %sw images. Fallback to defaults.\n',mpm_params.input(ccon).tag);
        MPMacq = hmri_get_defaults('MPMacq');
        mpm_params.input(ccon).TE = MPMacq.(['TE_' lower(mpm_params.input(ccon).tag) 'w']);
        mpm_params.input(ccon).TR = MPMacq.(['TR_' lower(mpm_params.input(ccon).tag) 'w']);
        mpm_params.input(ccon).fa = MPMacq.(['fa_' lower(mpm_params.input(ccon).tag) 'w']);
    end
end
  
% check that echo times are identical (common echoes only)
% NOTE: only necessary when not using the TE=0 extrapolation
if ~mpm_params.fullOLS
    for ccon = 1:mpm_params.ncon-1
        TEcon1 = mpm_params.input(ccon).TE;
        TEcon2 = mpm_params.input(ccon+1).TE;
        for necho = 1:min(length(TEcon1),length(TEcon2))
            if ~(TEcon1(necho)==TEcon2(necho))
                error('Echo times do not match between contrasts! Aborting.');
            end
        end
    end
end

% maximum TE for averaging (ms)
maxTEval4avg = 30; 
% find maximum number of echoes that are common to all available contrasts
ncommonTEvals = 1000;
for ccon = 1:mpm_params.ncon
    ncommonTEvals = min(length(mpm_params.input(ccon).TE),ncommonTEvals);
end
% find maximum number of echoes that are common to all available contrasts
% AND one more than the maxTEval4avg:
mpm_params.nr_echoes4avg = min(length(find(mpm_params.input(1).TE<maxTEval4avg))+1,ncommonTEvals);
fprintf(1,'\nINFO: averaged PDw/T1w/MTw will be calculated based on the first %d echoes.\n',mpm_params.nr_echoes4avg);
        
% if T1w and PDw data available, identify the protocol to define RF
% spoiling correction parameters (for T1 map calculation)
if mpm_params.PDidx && mpm_params.T1idx
    % retrieve all available protocols:
    MPMacq_sets = hmri_get_defaults('MPMacq_set');
    % current protocol is defined by [TR_pdw TR_t1w fa_pdw fa_t1w]:
    MPMacq_prot = [mpm_params.input(mpm_params.PDidx).TR;
                   mpm_params.input(mpm_params.T1idx).TR;
                   mpm_params.input(mpm_params.PDidx).fa;
                   mpm_params.input(mpm_params.T1idx).fa]';
    % then match the values and find protocol tag
    nsets = numel(MPMacq_sets.vals);
    ii = 0; mtch = false;
    while ~mtch && ii < nsets
        ii = ii+1;
        if all(MPMacq_prot == MPMacq_sets.vals{ii})
            mtch  = true;
            prot_tag = MPMacq_sets.tags{ii};
            fprintf(1,'\nINFO: MPM acquisition protocol = %s.\n', prot_tag);
        end
    end
    if ~mtch
        prot_tag = 'Unknown';
        fprintf(1,'\nWARNING: MPM protocol unknown. No RF spoiling correction will be applied.\n');
    end
else
    prot_tag = 'Unknown';
end
% now retrieve RF spoiling correction parameters
mpm_params.proc.RFC = hmri_get_defaults(['rfcorr.',prot_tag]);

% other processing parameters from defaults
% load threshold to save qMRI maps
mpm_params.proc.threshall = hmri_get_defaults('qMRI_maps_thresh');
% load PD maps processing parameters
mpm_params.proc.PD = hmri_get_defaults('PDproc');
if ~mpm_params.T1idx && mpm_params.proc.PD.PDmap
    fprintf(1,'\nWARNING: PD map calculation enabled but T1w images not available.\n\tPD map won''t be calculated.\n');
    mpm_params.proc.PD.PDmap = 0;
end
% whether OLS R2* is calculated
mpm_params.proc.R2sOLS = hmri_get_defaults('R2sOLS');

% consistency check for number of echoes averaged for A calculation:
if mpm_params.PDidx && mpm_params.T1idx 
    if mpm_params.proc.PD.nr_echoes_forA > size(mpm_params.input(mpm_params.T1idx).fnam,1)
        fprintf(1,['\nWARNING: number of T1w echoes to be averaged for PD calculation (%d)' ...
            '\nis bigger than the available number of echoes (%d). Setting nr_echoes_forA' ...
            '\nto the maximum number of echoes.\n'],mpm_params.proc.PD.nr_echoes_forA, ...
            size(mpm_params.input(mpm_params.T1idx).fnam,1));
        mpm_params.proc.PD.nr_echoes_forA = size(mpm_params.input(T1idx).fnam,1);
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