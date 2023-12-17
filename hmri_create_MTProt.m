function [fR1, fR2s, fMT, fA, PPDw, PT1w, PMTw]  = hmri_create_MTProt(jobsubj) %#ok<*STOUT>
%==========================================================================
% This is hmri_create_MTProt, part of the hMRI-Toolbox.
%
% PURPOSE
% Estimation of quantitative parameters (R1, R2*, PD, MT) from
% multi-contrast multi-echo FLASH protocol 
% 
% FORMAT
% [fR1, fR2s, fMT, fA, PPDw, PT1w, PMTw]  = hmri_create_MTProt(jobsubj)
%
% INPUTS
%   jobsubj     parameters for one subject out of the job list.
%               NB: ONE SINGLE DATA SET FROM ONE SINGLE SUBJECT IS
%               PROCESSED HERE, LOOP OVER SUBJECTS DONE AT HIGHER LEVEL.
%   P_trans     filenames of a pair of magnitude image + transmit bias map
%               [p.u.] of B1+ field (not mandatory), saved as
%               "b1_trans_input" field in jobsubj.
%
% OUTPUTS
%   fR1     R1 map output filename (only quantitative if B1+ map provided)
%   fR2s    R2* map output filename (simple fit on PD data or OLS across
%           all contrasts, according to defaults settings)
%   fMT     Magnetization transfer (MT) map output filename  
%   fA      Proton density map output filename (free water concentration
%           (PD) [%] or signal amplitude (A) [a.u.] depending on job and
%           defaults settings (PDproc & RF sensitivity bias correction)). 
%   PPDw    average PD-weighted image filename (or OLS fit at TE=0 if fullOLS = true) 
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
%   Imperfect spoiling correction
%       Preibisch and Deichmann, MRM 61:125-135 (2009)
%       Corrects for residual transverse magnetization coherences that
%       remain despite RF & gradient spoiling and result in deviations from
%       the Ernst equation (ideal expected signal in a FLASH acquisition). 
%       Modelling data and P2_a, P2_b correction factors calculation based
%       on code supplied by Ralf Deichmann. Only a small subset of
%       experimental parameters have been used and imperfect spoiling
%       correction can only be applied if these correction factors are
%       defined specifically for the acquisition protocol used.
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

flags = jobsubj.log.flags;
flags.PopUp = false;
hmri_log(sprintf('\t============ MPM PROCESSING - %s.m (%s) ============', mfilename, datetime('now')),flags);

% retrieves all required parameters for MPM processing
mpm_params = get_mpm_params(jobsubj);

% for convenience, define a few parameters to make formulae more readable
% and avoid number of repetitions:

% B1 transmit input:
% P_trans(1,:) = magnitude image (anatomical reference for coregistration)
% P_trans(2,:) = B1 map (p.u.)
P_trans = jobsubj.b1_trans_input;

% index number for each contrast - zero index means no images available
PDwidx = mpm_params.PDwidx;
T1widx = mpm_params.T1widx;
MTwidx = mpm_params.MTwidx;
% TE/TR/FA for each contrast 
% T1w and MTw contrasts are optional.
% PDw must be present or nothing can be calculated. The script will abort
% at the next line if no PDw echoes available. In that case, a warning has
% been thrown earlier (in get_mpm_params).
TE_pdw = mpm_params.input(PDwidx).TE;
TR_pdw = mpm_params.input(PDwidx).TR;
fa_pdw = mpm_params.input(PDwidx).fa;
if T1widx
    TE_t1w = mpm_params.input(T1widx).TE; %#ok<NASGU>
    TR_t1w = mpm_params.input(T1widx).TR;
    fa_t1w = mpm_params.input(T1widx).fa;
end
if MTwidx
    TE_mtw = mpm_params.input(MTwidx).TE; %#ok<NASGU>
    TR_mtw = mpm_params.input(MTwidx).TR;
    fa_mtw = mpm_params.input(MTwidx).fa;
end
% other parameters
threshall = mpm_params.proc.threshall;
PDproc = mpm_params.proc.PD;
ISC = mpm_params.proc.ISC;
dt = [spm_type('float32'),spm_platform('bigend')]; % for nifti output
outbasename = spm_file(mpm_params.input(PDwidx).fnam(1,:),'basename'); % for all output files
calcpath = mpm_params.calcpath;
mpm_params.outbasename = outbasename;
respath = mpm_params.respath;
supplpath = mpm_params.supplpath;
% Number of echoes averaged (maximum number or echoes available for ALL
% contrasts AND under TE_limit (+1) - see get_mpm_params)
avg_nr = mpm_params.nr_echoes4avg; 

% load PDw images
V_pdw = spm_vol(mpm_params.input(PDwidx).fnam);

% We set the PDw images as reference space for all results/different contrasts. 
% Therefore, the matrix dimension defined below (dm) is used across the whole script. 
% It should not be redefined.
V_ref = V_pdw(1);
dm = V_ref.dim;

%% =======================================================================%
% Calculate R2* map from PDw echoes
%=========================================================================%
fR2s = '';
hmri_log(sprintf('\t-------- R2* map calculation (basic exponential decay) --------'),mpm_params.nopuflags);
hmri_log(sprintf(['INFO: minimum number of echoes required for R2* map calculation is %d.' ...
    '\nNumber of PDw echoes available: %d'], mpm_params.neco4R2sfit, length(mpm_params.input(PDwidx).TE)), mpm_params.nopuflags);
if mpm_params.basicR2s
    if length(mpm_params.input(PDwidx).TE)<6
        hmri_log(sprintf(['GENERAL WARNING: deriving R2* map from an echo train including ' ...
            '\nfewer than 6 echoes has not been validated nor investigated. ' ...
            '\nFor a robust estimation, the minimum number of echoes required ' ...
            '\ndepends on many factors, amongst which: ' ...
            '\n\t- SNR/resolution' ...
            '\n\t- distribution/spacing between TEs: note that early echoes are more' ...
            '\n\t\taffected by the PDw contrast.' ...
            '\nInterpret results carefully, with in mind a possible lack of robustness ' ...
            '\nand reliability in the R2* estimation.']),mpm_params.defflags);
    end
        
    spm_progress_bar('Init',dm(3),'R2* fit','planes completed');
    
    % create nifti object for output R2* map
    fR2s = fullfile(calcpath,[outbasename '_' mpm_params.output(mpm_params.qR2s).suffix '.nii']);
    Ni = hmri_create_nifti(fR2s, V_ref, dt, mpm_params.output(mpm_params.qR2s).descrip{1});
       
    % Fit R2* decay
    for p = 1:dm(3)
        data = zeros([dm(1:2),numel(V_pdw)]);
        for cecho = 1:numel(V_pdw)
            % Allows for reslicing across TE
            data(:,:,cecho) = hmri_read_vols(V_pdw(cecho),V_ref,p,mpm_params.interp);
        end
        % Don't want log of < 1 => assuming typical dynamic range of Dicom
        % data, e.g. 12bit @TODO
        data = max(data, 1);
        R2s = hmri_calc_R2s(struct('data',data,'TE',TE_pdw),mpm_params.R2s_fit_method);
        
        Ni.dat(:,:,p) = max(min(R2s,threshall.R2s),-threshall.R2s)*1000; % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');
    
    % Set and write metadata
    input_files = mpm_params.input(PDwidx).fnam;
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = 'Basic R2* map';
    Output_hdr.history.output.units = 's-1';
    set_metadata(fR2s,Output_hdr,mpm_params.json);
else 
    hmri_log(sprintf('INFO: No (basic) R2* map will be calculated.\nInsufficient number of PDw echoes.'), mpm_params.defflags);
end

%% =======================================================================%
% Reading and averaging the images 
%=========================================================================%
hmri_log(sprintf('\t-------- Reading and averaging the images --------'),mpm_params.nopuflags);

% Average only first few echoes for increased SNR and fit T2* 
Pavg = cell(1,mpm_params.ncon);
for ccon=1:mpm_params.ncon % loop over available contrasts
    Pavg{ccon}  = fullfile(calcpath,[outbasename '_' mpm_params.input(ccon).tag 'w.nii']);
    V           = spm_vol(mpm_params.input(ccon).fnam);
    Ni = hmri_create_nifti(Pavg{ccon}, V(1), dt, ...
        sprintf('Averaged %sw images - %d echoes', mpm_params.input(ccon).tag, avg_nr));

    spm_progress_bar('Init',V(1).dim(3),Ni.descrip,'planes completed');
    for p = 1:V(1).dim(3)
        Y = zeros(V(1).dim(1:2));
        for nr = 1:avg_nr
            Y  = Y + hmri_read_vols(V(nr),V(1),p,mpm_params.interp);
        end
        Ni.dat(:,:,p) = Y/avg_nr;
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');
    
    input_files = mpm_params.input(ccon).fnam;
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = Ni.descrip;
    Output_hdr.history.output.units = 'a.u.';
    set_metadata(Pavg{ccon},Output_hdr,mpm_params.json);
end

% Average T1w image for PD calculation 
% (average over PDproc.nr_echoes_forA echoes, see hmri_defaults):
% @TODO - ensure reference volume is always consistent, e.g. V v's V_pdw to
% ensure we are always in the same reference space of the PDw volume.
if (PDwidx && T1widx)
    PT1w_forA = fullfile(calcpath,[outbasename '_T1w_forA.nii']);
    V           = spm_vol(mpm_params.input(T1widx).fnam);
    Ni = hmri_create_nifti(PT1w_forA, V(1), dt, ...
        sprintf('Averaged T1w images for PD calculation - %d echoes',PDproc.nr_echoes_forA));
 
    spm_progress_bar('Init',V(1).dim(3),Ni.descrip,'planes completed');
    for p = 1:V(1).dim(3)
        Y = zeros(V(1).dim(1:2));
        for nr = 1:PDproc.nr_echoes_forA
            Y  = Y + hmri_read_vols(V(nr),V(1),p,mpm_params.interp);
        end
        Ni.dat(:,:,p) = Y./PDproc.nr_echoes_forA;
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear');
end

%% =======================================================================%
% Coregistering the images
%=========================================================================%
hmri_log(sprintf('\t-------- Coregistering the images --------'),mpm_params.nopuflags);
% NOTE: coregistration can be disabled using the hmri_def.coreg2PDw flag

contrastCoregParams = zeros(mpm_params.ncon, 6);
if mpm_params.coreg
    if MTwidx
        contrastCoregParams(MTwidx,:) = hmri_coreg(Pavg{PDwidx}, Pavg{MTwidx}, mpm_params.coreg_flags);
    end
    if T1widx 
        contrastCoregParams(T1widx,:) = hmri_coreg(Pavg{PDwidx}, Pavg{T1widx}, mpm_params.coreg_flags);
        hmri_coreg(Pavg{PDwidx}, PT1w_forA, mpm_params.coreg_flags);
    end
end

V_trans = [];
if ~isempty(P_trans)
    % Load B1 mapping data if available and coregister to PDw
    % P_trans(1,:) = magnitude image (anatomical reference for coregistration) 
    % P_trans(2,:) = B1 map (p.u.)
    if mpm_params.coreg
        hmri_coreg(Pavg{PDwidx}, P_trans, mpm_params.coreg_bias_flags);
    end
    V_trans = spm_vol(P_trans);
end

% parameters saved for quality assessment
if mpm_params.QA.enable
    if exist(mpm_params.QA.fnam,'file')
        mpm_params.QA = spm_jsonread(mpm_params.QA.fnam);
    end
    if MTwidx
        mpm_params.QA.ContrastCoreg.MT2PD = contrastCoregParams(MTwidx,:);
    end
    if T1widx
        mpm_params.QA.ContrastCoreg.T12PD = contrastCoregParams(T1widx,:);
    end
    spm_jsonwrite(mpm_params.QA.fnam, mpm_params.QA, struct('indent','\t'));
end

% load averaged images
Vavg(PDwidx) = spm_vol(Pavg{PDwidx});
if MTwidx; Vavg(MTwidx) = spm_vol(Pavg{MTwidx}); end
if T1widx
    Vavg(T1widx) = spm_vol(Pavg{T1widx}); 
    VT1w_forA = spm_vol(PT1w_forA);
end


%% =======================================================================%
% multi-contrast R2* map calculation for quality assessment by AL
% METHOD: R2s map calculated for each contrast separately (PDw, MTw, T1w)
% to evaluate the SD of each R2* map in white matter (i.e. intra-run motion
% measurement).
%=========================================================================%
if mpm_params.QA.enable
    hmri_log(sprintf('\t-------- Multi-contrast R2* map calculation for QA --------'),mpm_params.nopuflags);
    
    fR2sQA = cell(1,mpm_params.ncon);
    for ccon = 1:mpm_params.ncon
        % only if enough echoes
        if mpm_params.estaticsR2s(ccon)
            
            fR2sQA{ccon} = fullfile(calcpath,[outbasename '_R2s_' mpm_params.input(ccon).tag 'w.nii']);
            Ni = hmri_create_nifti(fR2sQA{ccon}, V_ref, dt, ...
                'OLS R2* map [s-1]');
            
            Vcon = spm_vol(mpm_params.input(ccon).fnam);
            
            % The assumption is that the result of co-registering the average
            % weighted volumes is applicable for each of the echoes of that
            % contrast => Replicate the mat field across contrasts for all echoes.
            % % matField = cat(3, repmat(VPDw.mat, [1, 1, nPD]), ...
            % % repmat(VMTw.mat, [1, 1, nMT]), repmat(VT1w.mat, [1, 1, nT1]));
            spm_progress_bar('Init',dm(3),'multi-contrast R2* fit','planes completed');
            for p = 1:dm(3)
                data = zeros([dm(1:2),numel(Vcon)]);
                for cecho = 1:numel(mpm_params.input(ccon).TE)
                    % Allows for reslicing across TE
                    data(:,:,cecho) = hmri_read_vols(Vcon(cecho),V_ref,p,mpm_params.interp);
                end
                % Don't want log of < 1 => assuming typical dynamic range of Dicom
                % data, e.g. 12bit @TODO
                data = max(data, 1);
                R2s = hmri_calc_R2s(struct('data',data,'TE',mpm_params.input(ccon).TE),'ols');
                
                Ni.dat(:,:,p) = max(min(R2s,threshall.R2s),-threshall.R2s)*1000; % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), negative values are allowed to preserve Gaussian distribution
                spm_progress_bar('Set',p);
            end
            spm_progress_bar('Clear');
        end
    end
end

%% =======================================================================%
% OLS R2* map calculation by MFC
% [Reference: ESTATICS, Weiskopf et al. 2014]
%=========================================================================%
fR2s_OLS = '';
hmri_log(sprintf('\t-------- %s R2* map calculation (ESTATICS) --------',mpm_params.R2s_fit_method),mpm_params.nopuflags);
LogMsg = sprintf(['INFO: minimum number of echoes required for R2* map calculation is %d.' ...
    '\nNumber of echoes available: '], mpm_params.neco4R2sfit);
for ccon = 1:mpm_params.ncon
    LogMsg = sprintf('%s %sw = %d.',LogMsg, mpm_params.input(ccon).tag, length(mpm_params.input(ccon).TE));
end
hmri_log(LogMsg, mpm_params.nopuflags);

if mpm_params.proc.R2sOLS && any(mpm_params.estaticsR2s)
    if length(mpm_params.input(PDwidx).TE)<6 
        hmri_log(sprintf(['GENERAL WARNING: deriving R2* map from echo trains including ' ...
            '\nfewer than 6 echoes has not been validated nor investigated. ' ...
            '\nFor a robust estimation, the minimum number of echoes required ' ...
            '\ndepends on many factors, amongst which: ' ...
            '\n\t- SNR/resolution' ...
            '\n\t- distribution/spacing between TEs: note that early echoes are more' ...
            '\n\t\taffected by the specific contrast, violating the ESTATICS assumption ' ...
            '\n\t\tof a common decay between contrasts.' ...
            '\n\t- number of contrasts available (fewer echoes per contrast required ' ...
            '\n\t\tfor 3 (PDw, T1w, MTw) contrasts as compared to 2 or even 1) ' ...
            '\nInterpret results carefully, with in mind a possible lack of robustness ' ...
            '\nand reliability in the R2* estimation.']),mpm_params.defflags);
    end
        
    % OLS fit at TE=0: to be used instead of averaged echoes of each
    % contrast if "fullOLS" option is enabled
    if mpm_params.fullOLS
        hmri_log(sprintf('\t-------- and fit to TE=0 for each contrast --------'),mpm_params.nopuflags);
 
        Nmap = nifti;
        Pte0 = cell(1,mpm_params.ncon);
        for ccon = 1:mpm_params.ncon
            % requires a minimum of neco4R2sfit echoes for a robust fit
            if mpm_params.estaticsR2s(ccon)
                Pte0{ccon}  = fullfile(calcpath,[outbasename '_' mpm_params.input(ccon).tag 'w_' mpm_params.R2s_fit_method 'fit_TEzero.nii']);
                Ni = hmri_create_nifti(Pte0{ccon}, V_ref, dt, ...
                    sprintf('%s fit to TE=0 for %sw images - %d echoes', mpm_params.R2s_fit_method, mpm_params.input(ccon).tag, length(mpm_params.input(ccon).TE)));
                                
                % store processing history in metadata
                input_files = mpm_params.input(ccon).fnam;
                Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
                Output_hdr.history.output.imtype = Ni.descrip;
                Output_hdr.history.output.units = 'a.u.';

                % copy acquisition metadata so extrapolated data could be fed back into 
                % the toolbox if needed
                Output_hdr.acqpar = struct('RepetitionTime',mpm_params.input(ccon).TR, ...
                    'EchoTime',0, ...
                    'FlipAngle',mpm_params.input(ccon).fa);

                set_metadata(Pte0{ccon},Output_hdr,mpm_params.json);
                
                % re-load the updated NIFTI file (in case extended header has
                % been added, the offset has changed and must be updated before
                % writing the data to the file!)
                Nmap(ccon) = nifti(Pte0{ccon});
            else
                Pte0{ccon} = '';
            end
        end
    end % init nifti objects for fullOLS case
    
    fR2s_OLS    = fullfile(calcpath,[outbasename '_' mpm_params.output(mpm_params.qR2s).suffix '_' mpm_params.R2s_fit_method '.nii']);
    Ni = hmri_create_nifti(fR2s_OLS, V_ref, dt, ...
        [mpm_params.R2s_fit_method ' R2* map [s-1]']);
    
    % Combine the data and echo times:
    nechoes = zeros(1,mpm_params.ncon);
    for ccon = 1:mpm_params.ncon
        % only if enough echoes
        if mpm_params.estaticsR2s(ccon)
            nechoes(ccon) = size(mpm_params.input(ccon).fnam,1);
        end
    end
        
    spm_progress_bar('Init',dm(3),[mpm_params.R2s_fit_method ' R2* fit','planes completed']);
    for p = 1:dm(3)
        
        for ccon = 1:mpm_params.ncon
            
            % only if enough echoes
            if mpm_params.estaticsR2s(ccon)
                data = zeros([dm(1:2), nechoes(ccon)]);
                Vcon = spm_vol(mpm_params.input(ccon).fnam);
                
                for cecho = 1:nechoes(ccon)
                    % Additionally pass the result of co-registration since
                    % want output (e.g. intercept) in Vavg space.
                    % Don't want log of < 1 => assuming typical dynamic range of Dicom
                    % data, e.g. 12bit @TODO
                    data(:,:,cecho) = max(hmri_read_vols(Vcon(cecho),V_ref,p,mpm_params.interp, contrastCoregParams(ccon,:)), 1);
                    
                end
                dataToFit(ccon).data = data;
                dataToFit(ccon).TE = mpm_params.input(ccon).TE;
            end   
        end
        [R2s, intercepts] = hmri_calc_R2s(dataToFit,mpm_params.R2s_fit_method);

        % Writes "fullOLS" images (OLS fit to TE=0 for each contrast)
        if mpm_params.fullOLS
            for ccon = 1:mpm_params.ncon
                % only if enough echoes
                if mpm_params.estaticsR2s(ccon)
                    Nmap(ccon).dat(:,:,p) = intercepts{ccon};
                end
            end
        end
                
        % NB: mat field defined by V_pdw => first PDw echo
        % threshold T2* at +/- 0.1ms or R2* at +/- 10000 *(1/sec), 
        % negative values are allowed to preserve Gaussian distribution.
        Ni.dat(:,:,p) = max(min(R2s,threshall.R2s),-threshall.R2s)*1000; 
        spm_progress_bar('Set',p);
        
    end
    spm_progress_bar('Clear');

    % Set metadata (R2S_OLS)
    input_files = mpm_params.input(PDwidx).fnam;
    if (T1widx); input_files = char(input_files, mpm_params.input(T1widx).fnam); end
    if (MTwidx); input_files = char(input_files, mpm_params.input(MTwidx).fnam); end
    Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
    Output_hdr.history.output.imtype = mpm_params.output(mpm_params.qR2s).descrip;
    Output_hdr.history.output.units = mpm_params.output(mpm_params.qR2s).units;
    set_metadata(fR2s_OLS,Output_hdr,mpm_params.json);
   
else
    hmri_log(sprintf('INFO: No R2* map will be calculated using the ESTATICS model. \nInsufficient number of echoes.'), mpm_params.defflags);
end % OLS code

% if "fullOLS" option enabled AND could be applied to every contrast
% (enough echoes for every contrast), Pte0 images replace Pavg images in
% the rest of the code. We need to apply the substitution and reload the
% images... 
if mpm_params.fullOLS && all(mpm_params.estaticsR2s)
    for ccon = 1:mpm_params.ncon
        Pavg{ccon} = Pte0{ccon};
        Vavg(ccon) = spm_vol(Pavg{ccon});
    end
end
        

%% =======================================================================%
% Prepare output for R1, PD and MT maps
%=========================================================================%

Nmap    = nifti;
for ii=1:length(mpm_params.output)-1*(~isempty(fR2s)||~isempty(fR2s_OLS)) % R2s output already done
    
    % define NIFTI objects for output images
    outputFn = fullfile(calcpath,[outbasename '_' mpm_params.output(ii).suffix '.nii']);
    
    Nmap(ii) = hmri_create_nifti(outputFn, V_ref, dt, ...
        fullfile(calcpath,[outbasename '_' mpm_params.output(ii).suffix '.nii']));
    
end

fR1 = '';
fA = '';
fMT = '';
if (PDwidx && T1widx)
    fR1 = fullfile(calcpath,[outbasename '_' mpm_params.output(mpm_params.qR1).suffix '.nii']);
    fA  = fullfile(calcpath,[outbasename '_' mpm_params.output(mpm_params.qPD).suffix '.nii']);
    if MTwidx
        fMT = fullfile(calcpath,[outbasename '_'  mpm_params.output(mpm_params.qMT).suffix '.nii']);
    end        
end


%% =======================================================================%
% Map calculation continued (R1, PD, MT) 
%=========================================================================%
hmri_log(sprintf('\t-------- Map calculation continued (R1, MTR) --------'), mpm_params.nopuflags);

fa_pdw_rad = fa_pdw * pi / 180;
if MTwidx; fa_mtw_rad = fa_mtw * pi / 180; end
if T1widx; fa_t1w_rad = fa_t1w * pi / 180; end

spm_progress_bar('Init',dm(3),'Calculating maps','planes completed');

% First calculate R1 & MTR
for p = 1:dm(3)

    % PDw images are always available, so this bit is always loaded:
    PDw = hmri_read_vols(Vavg(PDwidx),Ni,p,mpm_params.interp);
    
    if ~isempty(V_trans)
        f_T = hmri_read_vols(V_trans(2,:),Ni,p,mpm_params.interp)/100; % divide by 100, since p.u. maps
    else
        f_T = [];
    end
    
    % Standard magnetization transfer ratio (MTR) in percent units [p.u.]
    % only if trpd = trmt and fapd = famt and if PDw and MTw images are
    % available
    if (MTwidx && PDwidx)
        MTw = hmri_read_vols(Vavg(MTwidx),Ni,p,mpm_params.interp);
        if (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw) % additional MTR image...
            MTR = (PDw-MTw)./(PDw+eps) * 100;
            % write MTR image
            Nmap(mpm_params.qMTR).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
        end          
    end
    
    % T1 map and A/PD maps can only be calculated if T1w images are
    % available:
    if T1widx

        T1w = hmri_read_vols(Vavg(T1widx),Ni,p,mpm_params.interp);
        
        % Transmit bias corrected quantitative T1 values; if f_T empty then semi-quantitative
        R1=hmri_calc_R1(struct('data',PDw,'fa',fa_pdw_rad,'TR',TR_pdw,'B1',f_T),...
            struct('data',T1w,'fa',fa_t1w_rad,'TR',TR_t1w,'B1',f_T),...
            mpm_params.small_angle_approx);
            
        if ISC.enabled&&~isempty(f_T)
            % MFC: We do have P2_a and P2_b parameters for this sequence
            % => T1 = A(B1) + B(B1)*T1app (see Preibisch 2009)
            R1 = R1./(...
                +(ISC.P2_a(1)*f_T.^2 + ISC.P2_a(2)*f_T + ISC.P2_a(3)).*R1 ...
                +(ISC.P2_b(1)*f_T.^2 + ISC.P2_b(2)*f_T + ISC.P2_b(3)) ...
                );
        end
        
        R1 = R1*1e6;
        R1(R1<0) = 0;
        tmp      = R1;
        tmp(isnan(tmp)) = 0;
        Nmap(mpm_params.qR1).dat(:,:,p) = min(max(tmp,-threshall.R1),threshall.R1)*1e-3; % truncating images
                
    end
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');


% apply UNICORT if required:
% NOTE: the (unicort) B1 map will then be used to calculate PD, but
% not to correct R1 map from RF spoiling effects (unsuitable
% interfering and higher order corrections!) nor to correct MT map
% from higher order transmit bias effects.
V_trans_unicort = [];
if (mpm_params.UNICORT.R1)
    out_unicort = hmri_create_unicort(Pavg{PDwidx}, fR1, jobsubj);
    fR1 = out_unicort.R1u{1};
    P_trans_unicort = out_unicort.B1u{1};
    V_trans_unicort = spm_vol(P_trans_unicort);
    Output_hdr = get_metadata(fR1);
    Output_hdr{1}.history.output.imtype = mpm_params.output(mpm_params.qR1).descrip;
    set_metadata(fR1,Output_hdr{1},mpm_params.json);
end

hmri_log(sprintf('\t-------- Map calculation continued (MT, A/PD) --------'), mpm_params.nopuflags);

spm_progress_bar('Init',dm(3),'Calculating maps (continued)','planes completed');

% Then calculate MT & PD
for p = 1:dm(3)

    if ~isempty(V_trans)
        f_T = hmri_read_vols(V_trans(2,:),Ni,p,mpm_params.interp)/100; % divide by 100, since p.u. maps
    elseif ~isempty(V_trans_unicort)
        f_T = hmri_read_vols(V_trans_unicort(1,:),Ni,p,mpm_params.interp)/100; % divide by 100, since p.u. maps        
    else
        f_T = [];
    end
    
    % PDw images are always available, so this bit is always loaded:
    PDw = hmri_read_vols(Vavg(PDwidx),Ni,p,mpm_params.interp);
    
    % T1 map and A/PD maps can only be calculated if T1w images are
    % available:
    if T1widx

        T1w = hmri_read_vols(Vavg(T1widx),Ni,p,mpm_params.interp);
        
        % if "fullOLS" option enabled, use the OLS fit at TE=0 as
        % "T1w_forA"; otherwise use the average calculated earlier (by
        % default, corresponds to the first echo to reduce R2* bias)
        if mpm_params.fullOLS
            T1w_forA = T1w;
        else
            T1w_forA = hmri_read_vols(VT1w_forA,Ni,p,mpm_params.interp);
        end
                
        % A values proportional to PD
        % f_T correction is applied either if:
        % - f_T has been provided as separate B1 mapping measurement (not
        % UNICORT!) or
        % - f_T has been calculated using UNICORT *AND* the UNICORT.PD flag
        % is enabled (advanced user only! method not validated yet!)
        if(~isempty(f_T))&&(~mpm_params.UNICORT.R1 || mpm_params.UNICORT.PD)
            A = hmri_calc_A(struct('data',PDw,'fa',fa_pdw_rad,'TR',TR_pdw,'B1',f_T),...
                struct('data',T1w_forA,'fa',fa_t1w_rad,'TR',TR_t1w,'B1',f_T),...
                mpm_params.small_angle_approx);
        else
            % semi-quantitative A
            A = hmri_calc_A(struct('data',PDw,'fa',fa_pdw_rad,'TR',TR_pdw,'B1',1),...
                struct('data',T1w_forA,'fa',fa_t1w_rad,'TR',TR_t1w,'B1',1),...
                mpm_params.small_angle_approx);
        end
        
        tmp      = A;
        tmp(isinf(tmp)) = 0;
        tmp(isnan(tmp)) = 0;
        Nmap(mpm_params.qPD).dat(:,:,p) = max(min(tmp,threshall.A),-threshall.A);
        % dynamic range increased to 10^5 to accommodate phased-array coils and symmetrical for noise distribution

        % for MT maps calculation, one needs MTw images on top of the T1w
        % and PDw ones...
        if MTwidx
            MTw = hmri_read_vols(Vavg(MTwidx),Ni,p,mpm_params.interp);

            switch mpm_params.MTsatB1CorrectionModel
                case 'helms' 
                    B1_mtw=1;
                    if isempty(f_T)
                        hmri_log('WARNING: MTsat B1-correction using the Helms model was selected but no B1 map data was found. MTsat will only be corrected with the quadratic model from Helms, et al. (MRM 2008).', mpm_params.defflags)
                    end
                case 'lipp'
                    B1_mtw=f_T;
                    if isempty(f_T)
                        hmri_log('WARNING: MTsat B1-correction using the Lipp model was selected but no B1 map data was found. MTsat will not be B1 corrected.', mpm_params.defflags)
                    end
            end

            % MT in [p.u.]
            A_forMT = hmri_calc_A(struct('data',PDw,'fa',fa_pdw_rad,'TR',TR_pdw,'B1',B1_mtw),...
                struct('data',T1w_forA,'fa',fa_t1w_rad,'TR',TR_t1w,'B1',B1_mtw),...
                mpm_params.small_angle_approx);
            R1_forMT = hmri_calc_R1(struct('data',PDw,'fa',fa_pdw_rad,'TR',TR_pdw,'B1',B1_mtw),...
                struct('data',T1w,'fa',fa_t1w_rad,'TR',TR_t1w,'B1',B1_mtw),...
                mpm_params.small_angle_approx);
            MT = hmri_calc_MTsat(struct('data',MTw,'fa',fa_mtw_rad,'TR',TR_mtw,'B1',B1_mtw), A_forMT, R1_forMT);

            % f_T correction is applied either if:
            % - f_T has been provided as separate B1 mapping measurement (not
            % UNICORT!) or
            % - f_T has been calculated using UNICORT *AND* the UNICORT.MT flag
            % is enabled (advanced user only! method not validated yet!)
            if (~isempty(f_T))&&(~mpm_params.UNICORT.R1 || mpm_params.UNICORT.MT)
                MT = hmri_correct_MTsat(MT,f_T,mpm_params.MTsatB1CorrectionModel,mpm_params.MTsatB1CorrectionC);
            end

            MT(isnan(MT))=0;
            Nmap(mpm_params.qMT).dat(:,:,p) = max(min(MT,threshall.MT),-threshall.MT);
        end
    
    end

    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');


% set metadata for all output images
input_files = mpm_params.input(PDwidx).fnam;
if (T1widx); input_files = char(input_files, mpm_params.input(T1widx).fnam); end
if (MTwidx); input_files = char(input_files, mpm_params.input(MTwidx).fnam); end
Output_hdr = init_mpm_output_metadata(input_files, mpm_params);
for ctr = 1:length(mpm_params.output)-1
    Output_hdr.history.output.imtype = mpm_params.output(ctr).descrip;
    Output_hdr.history.output.units = mpm_params.output(ctr).units;
    set_metadata(fullfile(calcpath,[outbasename '_' mpm_params.output(ctr).suffix '.nii']),Output_hdr,mpm_params.json);
end

%% =======================================================================%
% ACPC Realign all images - only if MT map created
%=========================================================================%
if mpm_params.ACPCrealign 
    if (MTwidx && PDwidx && T1widx)
    hmri_log(sprintf('\t-------- ACPC Realign all images --------'), mpm_params.nopuflags);
    
    % Define and calculate masked MT image
    % Load MT image
    V_MT = spm_vol(fMT);
    MTimage = spm_read_vols(V_MT);    
    % Define new file name for masked MT image
    V_MT.fname = fullfile(calcpath,['masked_' spm_str_manip(fMT,'t')]);
    % Load average PDw image (mask based on averaged PDw image)
    PDWimage = spm_read_vols(Vavg(PDwidx));
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
    end
    
    % Save transformation matrix
    spm_jsonwrite(fullfile(supplpath,'hMRI_map_creation_ACPCrealign_transformation_matrix.json'),R,struct('indent','\t'));

    else
        hmri_log(sprintf(['WARNING: ACPC Realign was enabled, but no MT map was available \n' ...
                   'to proceed. ACPC realignment must be done separately, e.g. you can \n'...
                   'run [hMRI tools > Auto-reorient] before calculating the maps.\n' ...
                   'NOTE: segmentation might crash if no initial reorientation.']), mpm_params.defflags);
    end
end

% for quality assessment and/or PD map calculation
% segmentation preferentially performed on MT map but can be done on R1 map
% if no MT map available. Therefore, we must at least have R1 available,
% i.e. both PDw and T1w inputs...
if (mpm_params.QA.enable||(PDproc.calibr)) && (PDwidx && T1widx)
    if ~isempty(fMT)
        Vsave = spm_vol(fMT);
        threshMT=threshall.MT;
    else % ~isempty(fR1); 
        Vsave = spm_vol(fR1); 
        threshMT=threshall.R1*1e3;
    end
    MTtemp = spm_read_vols(Vsave);
    
    % The 5 outer voxels in all directions are nulled in order to remove
    % artefactual effects from the MT map on segmentation: 
    MTtemp(1:5,:,:)=0; MTtemp(end-5:end,:,:)=0;
    MTtemp(:,1:5,:)=0; MTtemp(:,end-5:end,:)=0;
    MTtemp(:,:,1:5)=0; MTtemp(:,:,end-5:end)=0;
    
    % Null very bright and negative voxels
    MTtemp(abs(MTtemp)==threshMT)=0;
    MTtemp(isinf(MTtemp))=0;
    MTtemp(isnan(MTtemp))=0;
    MTtemp(MTtemp<0)=0;
    
    Vsave.fname = spm_file(Vsave.fname,'suffix','_outer_suppressed');
    spm_write_vol(Vsave,MTtemp);
    
    % use unified segmentation with uniform defaults across the toobox:
    job_brainmask = hmri_get_defaults('segment');
    job_brainmask.channel.vols = {Vsave.fname};
    job_brainmask.channel.write = [0 0]; % no need to write BiasField nor BiasCorrected image
    output_list = spm_preproc_run(job_brainmask);
    fTPM = char(cat(1,output_list.tiss.c));
end

% for quality assessment - the above segmentation must have run
if mpm_params.QA.enable && exist('fTPM','var') && any(mpm_params.estaticsR2s)
    
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
        R2s = spm_read_vols(spm_vol(fR2sQA{ccon}));
        MaskedR2s = squeeze(R2s.*WMmask);
        SDR2s = std(MaskedR2s(MaskedR2s~=0),[],1);
        mpm_params.QA.SDR2s.([mpm_params.input(ccon).tag 'w']) = SDR2s;
    end
    
    spm_jsonwrite(mpm_params.QA.fnam, mpm_params.QA, struct('indent','\t'));
end

% PD map calculation
% for correction of the R2s bias in the A map if that option is enabled...
if PDproc.T2scorr && (~isempty(fR2s)||~isempty(fR2s_OLS))
    % uses OLS it if available - less noisy
    if ~isempty(fR2s_OLS)
        PR2s = fR2s_OLS;
    else
        PR2s = fR2s;
    end
    
    % calculate correction (expected to be between 1 and 1.5 approx)
    R2s = spm_read_vols(spm_vol(PR2s));
    R2scorr4A = zeros(size(R2s));
    for cecho=1:mpm_params.proc.PD.nr_echoes_forA
        TE = mpm_params.input(PDwidx).TE(cecho)*0.001; % in seconds
        R2scorr4A = R2scorr4A + exp(-TE.*R2s);
    end
    R2scorr4A = R2scorr4A/mpm_params.proc.PD.nr_echoes_forA;
    
    % save correction for inspection
    fR2scorr4A = spm_file(PR2s,'suffix','_corr4A');   
    NiR2scorr4A = hmri_create_nifti(fR2scorr4A, V_ref, dt, ...
        'R2* bias correction factor for A map (T2scorr option)');
   
    NiR2scorr4A.dat(:,:,:) = R2scorr4A;
    
     
    % apply correction
    fAcorr = spm_file(fA,'suffix','_R2scorr');
    NiAcorr = hmri_create_nifti(fAcorr, V_ref, dt, ...
        'R2* bias corrected A map (T2scorr option)');
    
    tmp = spm_read_vols(spm_vol(fA))./(R2scorr4A+eps);
    tmp(isnan(tmp)|isinf(tmp)) = 0;
    tmp = max(min(tmp,threshall.A),-threshall.A);
    NiAcorr.dat(:,:,:) = tmp;
    fA = fAcorr;
end

% PD map calculation continued
if ~isempty(f_T) && ~isempty(fA) && exist('fTPM','var') && (mpm_params.UNICORT.PD || ~mpm_params.UNICORT.R1)
    % if calibration enabled, do the Unified Segmentation bias correction
    % if required and calibrate the PD map
    if PDproc.calibr
        PDcalculation(fA,fTPM,mpm_params);
    end
end

% copy final result files into Results directory
% NB: to avoid ambiguity for users, only the 4 final maps to be used in
% further analysis are in Results, all other files (additional maps, 
% processing parameters, etc) are in Results/Supplementary.
if ~isempty(fR1)
    fR1_final = fullfile(respath, spm_file(fR1,'filename'));
    try copyfile(fR1,fR1_final); end
    try copyfile([spm_str_manip(fR1,'r') '.json'],[spm_str_manip(fR1_final,'r') '.json']); end %#ok<*TRYNC>
    fR1 = fR1_final;
end

if ~isempty(fR2s)
    fR2s_final = fullfile(respath, spm_file(fR2s,'filename'));
    copyfile(fR2s,fR2s_final);
    try copyfile([spm_str_manip(fR2s,'r') '.json'],[spm_str_manip(fR2s_final,'r') '.json']); end
    fR2s = fR2s_final;
end

% NB: if OLS calculation of R2s map has been done, the output file for R2s
% map is the OLS result. In that case, the basic R2s map is moved to
% Results/Supplementary while the R2s_OLS is copied into results directory: 
if mpm_params.proc.R2sOLS && ~isempty(fR2s_OLS)
    % move basic R2s map to Results/Supplementary
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

PPDw_final = fullfile(supplpath, spm_file(Pavg{PDwidx},'filename'));
copyfile(Pavg{PDwidx},PPDw_final);
try copyfile([spm_str_manip(Pavg{PDwidx},'r') '.json'],[spm_str_manip(PPDw_final,'r') '.json']); end
PPDw = PPDw_final;

if T1widx
    PT1w_final = fullfile(supplpath, spm_file(Pavg{T1widx},'filename'));
    copyfile(Pavg{T1widx},PT1w_final);
    try copyfile([spm_str_manip(Pavg{T1widx},'r') '.json'],[spm_str_manip(PT1w_final,'r') '.json']); end
    PT1w = PT1w_final;
else
    PT1w = '';
end

if MTwidx
    PMTw_final = fullfile(supplpath, spm_file(Pavg{MTwidx},'filename'));
    copyfile(Pavg{MTwidx},PMTw_final);
    try copyfile([spm_str_manip(Pavg{MTwidx},'r') '.json'],[spm_str_manip(PMTw_final,'r') '.json']); end
    PMTw = PMTw_final;
else
    PMTw = '';
end

% save processing params (mpm_params)
spm_jsonwrite(fullfile(supplpath,'hMRI_map_creation_mpm_params.json'),mpm_params,struct('indent','\t'));

spm_progress_bar('Clear');

hmri_log(sprintf('\t============ MPM PROCESSING: completed (%s) ============', datetime('now')),mpm_params.nopuflags);


end


%% =======================================================================%
% Proton density map calculation
%=========================================================================%
function PDcalculation(fA, fTPM, mpm_params)
% fA is the filename of the output A image 
% (not yet properly quantitative PD map when it enters PDcalculation)
% fTMP is the list of TPMs generated from MT map
hmri_log(sprintf('\t-------- Proton density map calculation --------'), mpm_params.nopuflags);

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
maskedA(isinf(maskedA)) = 0;
maskedA(isnan(maskedA)) = 0;
maskedA(maskedA==threshA) = 0;
maskedA(maskedA<0) = 0;
spm_write_vol(V_maskedA,maskedA);

% Bias-field correction of masked A map
% use unified segmentation with uniform defaults across the toolbox:
job_bfcorr = hmri_get_defaults('segment');
job_bfcorr.channel.vols = {V_maskedA.fname};
job_bfcorr.channel.biasreg = PDproc.biasreg;
job_bfcorr.channel.biasfwhm = PDproc.biasfwhm;
job_bfcorr.channel.write = [1 0]; % need the BiasField, obviously!
for ctis=1:length(job_bfcorr.tissue)
    job_bfcorr.tissue(ctis).native = [0 0]; % no need to write c* volumes
end
output_list = spm_preproc_run(job_bfcorr);

% Bias field correction of A map.
% Bias field calculation is based on the masked A map, while correction
% must be applied to the unmasked A map. The BiasField is therefore
% retrieved from previous step and applied onto the original A map.
BFfnam = output_list.channel.biasfield{1};
BF = double(spm_read_vols(spm_vol(BFfnam)));
Y = BF.*spm_read_vols(spm_vol(fA));

% Calibration of flattened A map to % water content using typical white
% matter value from the hmri_defaults (see hmri_def.PDproc.WMval)
A_WM = WMmask.*Y; 
Y = Y/mean(A_WM(A_WM~=0))*PDproc.WMval;
hmri_log(sprintf(['INFO (PD calculation):\n\tmean White Matter intensity: %.1f\n' ...
    '\tSD White Matter intensity %.1f\n'],mean(A_WM(A_WM~=0)),std(A_WM(A_WM~=0))), mpm_params.defflags);
Y(Y>200) = 0;
% MFC: Estimating Error for data set to catch bias field issues:
errorEstimate = std(A_WM(A_WM > 0))./mean(A_WM(A_WM > 0));
Vsave = spm_vol(fA);
Vsave.descrip = [Vsave.descrip '. Error Estimate: ', num2str(errorEstimate)];
if errorEstimate > 0.06 %#ok<BDSCI>
    % MFC: Testing on 15 subjects showed 6% is a good cut-off:
    hmri_log(sprintf(['WARNING: Error estimate is high (%.1f%%) for calculated PD map:\n%s' ...
        '\nError higher than 6%% may indicate motion.\n'], errorEstimate*100, Vsave.fname), mpm_params.defflags);
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

% save SPM version (slight differences may appear in the results depending
% on the SPM version!)
[v,r] = spm('Ver');
mpm_params.SPMver = sprintf('%s (%s)', v, r);

% flags for logging information and warnings
mpm_params.defflags = jobsubj.log.flags; % default flags
mpm_params.nopuflags = jobsubj.log.flags; % force no Pop-Up
mpm_params.nopuflags.PopUp = false; 

hmri_log(sprintf('\t------------ PROCESSING PARAMETERS: SETTING UP AND CONSISTENCY CHECK ------------'),mpm_params.nopuflags);

% global parameters
mpm_params.json = hmri_get_defaults('json');
mpm_params.centre = hmri_get_defaults('centre');
mpm_params.calcpath = jobsubj.path.mpmpath;
mpm_params.respath = jobsubj.path.respath;
mpm_params.supplpath = jobsubj.path.supplpath;
mpm_params.QA.enable = hmri_get_defaults('qMRI_maps.QA'); % quality assurance
if mpm_params.QA.enable
    mpm_params.QA.fnam = fullfile(mpm_params.supplpath,'hMRI_map_creation_quality_assessment.json');
    spm_jsonwrite(mpm_params.QA.fnam,mpm_params.QA,struct('indent','\t'));
end
mpm_params.ACPCrealign = hmri_get_defaults('qMRI_maps.ACPCrealign'); % realigns qMRI maps to MNI
mpm_params.interp = hmri_get_defaults('interp');
mpm_params.fullOLS = hmri_get_defaults('fullOLS'); % uses all echoes for OLS fit at TE=0

mpm_params.neco4R2sfit = hmri_get_defaults('neco4R2sfit'); % minimum number of echoes for R2* calculation 
if mpm_params.neco4R2sfit<2
    hmri_log(sprintf(['WARNING: the (strict) minimum number of echoes required to ' ...
        '\ncalculate R2* map is 2. The default value (%d) has been modified ' ...
        '\nand is now neco4R2sfit = 2.']),mpm_params.defflags);
    mpm_params.neco4R2sfit = 2;
end

% UNICORT settings:
mpm_params.UNICORT.R1 = isfield(jobsubj.b1_type,'UNICORT'); % uses UNICORT to estimate B1 transmit bias
tmp = hmri_get_defaults('UNICORT');
mpm_params.UNICORT.PD = tmp.PD; % uses B1map estimated as biasfield for R1 to correct for B1 transmit bias in PD
mpm_params.UNICORT.MT = tmp.MT; % uses B1map estimated as biasfield for R1 to correct for B1 transmit bias in MT

LogMsg = '';
if mpm_params.UNICORT.R1
    LogMsg = sprintf(['%s\nINFO: B1 transmit field estimated using UNICORT ' ...
        '\nand applied to correct R1 maps for B1 transmit bias.'], LogMsg);
    flags = mpm_params.nopuflags;
end
if mpm_params.UNICORT.PD
    LogMsg = sprintf(['%s\nWARNING: UNICORT B1+ estimate (if available) will be ' ...
        '\nused to correct the PD map for B1 transmit bias. This method has ' ...
        '\nnot been validated. Use with care!'], LogMsg);
    flags = mpm_params.defflags;
end
if mpm_params.UNICORT.MT
    LogMsg = sprintf(['%s\nWARNING: UNICORT B1+ estimate (if available) will be ' ...
        '\nused to correct the MT map for high-order B1 transmit bias. This ' ...
        '\nmethod has not been validated. Use with care!'], LogMsg);
    flags = mpm_params.defflags;
end
if ~isempty(LogMsg)
    hmri_log(LogMsg,flags);
end

% check consistency in UNICORT settings:
if ~mpm_params.UNICORT.R1
    if mpm_params.UNICORT.PD
        hmri_log(sprintf(['WARNING: no UNICORT B1 estimate available for PD calculation.' ...
            '\nB1 bias correction must be defined as "UNICORT" (not the case).' ...
            '\nUNICORT.PD is disabled!']),mpm_params.defflags);
        mpm_params.UNICORT.PD = false;
    end
    if mpm_params.UNICORT.MT
        hmri_log(sprintf(['WARNING: no UNICORT B1 estimate available for MT calculation.' ...
            '\nB1 bias correction must be defined as "UNICORT" (not the case).' ...
            '\nUNICORT.MT is disabled!']),mpm_params.defflags);
        mpm_params.UNICORT.MT = false;
    end
end

% retrieve input file names for map creation.
% the "mpm_params.input" field is an array, each element corresponds to a
% contrast.  
% if no input files for a given contrast, no input entry created and
% warning is thrown.
ccon = 0;
LogMsg = 'INFO: FLASH echoes loaded for each contrast are: ';
% 1) try MTw contrast:
tmpfnam   = char(jobsubj.raw_mpm.MT); % P_mtw
if isempty(tmpfnam)
    LogMsg = sprintf('%s\n\t- WARNING: no MT-weighted FLASH echoes available!',LogMsg);
    mpm_params.MTwidx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    LogMsg = sprintf('%s\n\t- MT-weighted: %d echoes',LogMsg,size(tmpfnam,1));
    mpm_params.input(ccon).fnam = tmpfnam;
    mpm_params.input(ccon).tag = 'MT';  
    mpm_params.MTwidx = ccon;
end  
% 2) try PDw contrast:
tmpfnam   = char(jobsubj.raw_mpm.PD); % P_pdw
if isempty(tmpfnam)
    LogMsg = sprintf('%s\n\t- WARNING: no PD-weighted FLASH echoes available! \n\t\tThe map creation won''t proceed!',LogMsg);
    mpm_params.PDwidx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    LogMsg = sprintf('%s\n\t- PD-weighted: %d echoes',LogMsg,size(tmpfnam,1));
    mpm_params.input(ccon).fnam = tmpfnam;
    mpm_params.input(ccon).tag = 'PD';  
    mpm_params.PDwidx = ccon;
end  
% 3) try T1w contrast:
tmpfnam   = char(jobsubj.raw_mpm.T1); % P_t1w
if isempty(tmpfnam)
    LogMsg = sprintf('%s\n\t- WARNING: no T1-weighted FLASH echoes available!',LogMsg);
    mpm_params.T1widx = 0; % zero index means no contrast available
else
    ccon = ccon+1;
    LogMsg = sprintf('%s\n\t- T1-weighted: %d echoes',LogMsg,size(tmpfnam,1));
    mpm_params.input(ccon).fnam = tmpfnam;
    mpm_params.input(ccon).tag = 'T1';  
    mpm_params.T1widx = ccon; % zero index means no contrast available    
end 
mpm_params.ncon = ccon; % number of contrasts available

% Message displayed as pop-up if enabled since it is important information 
hmri_log(LogMsg, mpm_params.defflags);

% Give a warning if the contrasts PDw, MTw and T1w have non-matching
% dimensions or orientations
raw = [jobsubj.raw_mpm.MT; jobsubj.raw_mpm.PD; jobsubj.raw_mpm.T1];
[orientationsMatch, orientationsWarning] = hmri_check_nifti_orientations(spm_vol(char(raw)));
if ~orientationsMatch
    % Message displayed as pop-up if enabled since it is important information 
    LogMsg = sprintf(['WARNING: MTw, PDw and T1w images have non-matching dimensions', ...
        '\n\t- see console for more information']);
    hmri_log(LogMsg, mpm_params.defflags);

    % Message is verbose and so only printed to console and log
    hmri_log(orientationsWarning, setfield(mpm_params.defflags,'PopUp',false));
end

% collect TE, TR and FA for each available contrast
for ccon = 1:mpm_params.ncon
    p = get_trtefa(mpm_params.input(ccon).fnam);
    if ~all(cellfun('isempty',{p.tr p.te p.fa}))
        mpm_params.input(ccon).TE = cat(1,p.te);
        mpm_params.input(ccon).TR = p(1).tr;
        mpm_params.input(ccon).fa = p(1).fa;
    else
        hmri_log(sprintf('WARNING: No TE/TR/FA values found for %sw images. Fallback to defaults.',mpm_params.input(ccon).tag),mpm_params.defflags);
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
                hmri_log('ERROR: Echo times do not match between contrasts! Aborting.',mpm_params.defflags);
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
hmri_log(sprintf('INFO: averaged PDw/T1w/MTw will be calculated based on the first %d echoes.',mpm_params.nr_echoes4avg),mpm_params.nopuflags);

% parameter determining whether to use small flip angle assumption or not
mpm_params.small_angle_approx = hmri_get_defaults('small_angle_approx');
if mpm_params.PDwidx && mpm_params.T1widx 
    if mpm_params.small_angle_approx
        hmri_log('INFO: Calculation of R1 and PD will use the small angle approximation.',mpm_params.defflags);
    elseif ~mpm_params.small_angle_approx
        hmri_log('INFO: Calculation of R1 and PD will not use the small angle approximation.',mpm_params.defflags);
    else
        hmri_log(sprintf(['WARNING: mpm_params.small_angle_approx flag not set.' ...
            '\nPlease check that hmri_def.small_angle_approx is set' ...
            '\nin the appropriate defaults file.']),mpm_params.defflags);
    end
end

% if T1w and PDw data available, identify the protocol to define imperfect 
% spoiling correction parameters (for T1 map calculation)
ISC = hmri_get_defaults('imperfectSpoilCorr.enabled');
if ~ISC
    hmri_log(sprintf(['INFO: Imperfect spoiling correction is disabled.' ...
        '\nIf your data were acquired with one of the standard MPM ' ...
        '\nprotocols (customized MT-FLASH sequence) for which the correction ' ...
        '\ncoefficients are available, it is recommended to enable that option.']),mpm_params.defflags);
end
if mpm_params.PDwidx && mpm_params.T1widx && ISC
    % retrieve all available protocols:
    MPMacq_sets = hmri_get_defaults('MPMacq_set');
    % current protocol is defined by [TR_pdw TR_t1w fa_pdw fa_t1w]:
    MPMacq_prot = [mpm_params.input(mpm_params.PDwidx).TR;
                   mpm_params.input(mpm_params.T1widx).TR;
                   mpm_params.input(mpm_params.PDwidx).fa;
                   mpm_params.input(mpm_params.T1widx).fa]';
    % then match the values and find protocol tag
    nsets = numel(MPMacq_sets.vals);
    ii = 0; mtch = false;
    while ~mtch && ii < nsets
        ii = ii+1;
        % don't check exact equality, but equality to within certain amount
        % of floating point error
        if all(abs(MPMacq_prot - MPMacq_sets.vals{ii}) < 5*eps(MPMacq_prot))
            mtch  = true;
            prot_tag = MPMacq_sets.tags{ii};
            hmri_log(sprintf(['INFO: MPM acquisition protocol = %s.' ...
                '\n\tThe coefficients corresponding to this protocol will be applied' ...
                '\n\tto correct for imperfect spoiling. Please check carefully that' ...
                '\n\tthe protocol used is definitely the one for which the ' ...
                '\n\tcoefficients have been calculated.'],prot_tag),mpm_params.defflags);
        end
    end
    if ~mtch
        prot_tag = 'Unknown';
        hmri_log(sprintf(['WARNING: MPM protocol unknown. ' ...
            '\n\tCorrection for imperfect spoiling will not be applied.']),mpm_params.defflags);
        hmri_get_defaults('imperfectSpoilCorr.enabled',false);
        ISC = false;
    end
else
    prot_tag = 'Unknown';
end
% now retrieve imperfect spoiling correction coefficients
mpm_params.proc.ISC = hmri_get_defaults(['imperfectSpoilCorr.',prot_tag]);
if ISC % only check for small_angle_approx assumptions if ISC applied
    % assume ISC parameters without small_angle_approx flag were generated with small angle approximation
    if ~isfield(mpm_params.proc.ISC,'small_angle_approx')
        hmri_log(sprintf(['WARNING: Small angle approx flag not set for' ...
            '\nthe imperfect spoiling correction parameters' ...
            '\nassociated with this MPM protocol.' ...
            '\nAssuming they were generated using the small angle approximation.']),mpm_params.defflags);
        mpm_params.proc.ISC.small_angle_approx = true;
    end
    % check small angle approximation used consistently, as causes change in
    % correction parameters
    if (mpm_params.proc.ISC.small_angle_approx)&&(~mpm_params.small_angle_approx)
        hmri_log(sprintf(['WARNING: Imperfect spoiling correction (ISC) coefficients were computed using the' ...
            '\nsmall angle approximation but maps will be created without this approximation.' ...
            '\nCorrection for imperfect spoiling will thus not be applied.' ...
            '\nTo use ISC, either calculate ISC coefficients without the small angle' ...
            '\napproximation or set hmri_def.small_angle_approx = true in the appropriate' ...
            '\nlocal defaults file.']),mpm_params.defflags);
        mpm_params.proc.ISC.enabled = false;
    elseif (~mpm_params.proc.ISC.small_angle_approx)&&(mpm_params.small_angle_approx)
        hmri_log(sprintf(['WARNING: Imperfect spoiling correction (ISC) coefficients were computed without the' ...
            '\nsmall angle approximation but maps will be created using this approximation.' ...
            '\nCorrection for imperfect spoiling will thus not be applied.' ...
            '\nTo use ISC, either calculate ISC coefficients using the small angle' ...
            '\napproximation or set hmri_def.small_angle_approx = false in the appropriate' ...
            '\nlocal defaults file.']),mpm_params.defflags);
        mpm_params.proc.ISC.enabled = false;
    else
        mpm_params.proc.ISC.enabled = true;
    end
else
    mpm_params.proc.ISC.enabled = false;
end

% RF sensitivity bias correction
mpm_params.proc.RFsenscorr = jobsubj.sensitivity;

% other processing parameters from defaults
% load threshold to save qMRI maps
mpm_params.proc.threshall = hmri_get_defaults('qMRI_maps_thresh');
% load PD maps processing parameters (including calibr (calibration
% parameter) and T2scorr (T2s correction) fields)
mpm_params.proc.PD = hmri_get_defaults('PDproc'); %todo: check user supplied parameters.
% if no RF sensitivity bias correction or no B1 transmit bias correction
% applied, not worth trying any calibration:
if (isfield(mpm_params.proc.RFsenscorr,'RF_none')||(isempty(jobsubj.b1_trans_input)&&~mpm_params.UNICORT.PD)) && mpm_params.proc.PD.calibr
    hmri_log(sprintf(['WARNING: both RF sensitivity and B1 transmit bias corrections '...
        '\nare required to generate a quantitative (calibrated) PD map.' ...
        '\nEither or both of these is missing. Therefore an amplitude ' ...
        '\n"A" map will be output instead of a quantitative ' ...
        '\nPD map. PD map calibration has been disabled.']),mpm_params.defflags);
    mpm_params.proc.PD.calibr = 0;
end    

% if fullOLS, T2*-weighting bias correction must not be applied
if mpm_params.fullOLS && mpm_params.proc.PD.T2scorr
    hmri_log(sprintf(['WARNING: if TE=0 fit is enabled (fullOLS option), ' ...
        '\nno T2*-weighting bias correction is required. \nT2scorr option disabled.']),mpm_params.defflags);
    mpm_params.proc.PD.T2scorr = 0;
end   
% T2*-weighting bias correction must always be applied somehow. If neither
% the fullOLS nor the T2scorr options are enabled, we don't force it to be
% applied (could still be usefull for comparison purposes) but throw a
% warning:
if ~mpm_params.fullOLS && ~mpm_params.proc.PD.T2scorr
    hmri_log(sprintf(['WARNING: both TE=0 fit (fullOLS option) and T2*-weighting ' ...
        '\nbias correction (T2scorr option) are disabled. The T2* bias won''t ' ...
        '\nbe corrected for and impact the resulting PD map. It is strongly ' ...
        '\nrecommended to have either of these options enabled!' ...
        '\nNOTE: this recommendation does not apply to single-echo datasets.']),mpm_params.defflags);
end   
    
% whether OLS R2* is calculated
mpm_params.proc.R2sOLS = hmri_get_defaults('R2sOLS');

% method to calculate R2*
mpm_params.R2s_fit_method = hmri_get_defaults('R2s_fit_method');
if mpm_params.proc.R2sOLS
    outstring=sprintf('Using %s method to compute R2*.',mpm_params.R2s_fit_method);
    switch lower(mpm_params.R2s_fit_method)
        case {'ols'}
             outstring=[outstring sprintf('\nNote that %s is sub-optimal, and that "wls1" might improve your R2* maps at the expense of a longer computation (Edwards et al. ISMRM 2022).',mpm_params.R2s_fit_method)];
        case {'wls1','wls2','wls3'}
             outstring=[outstring sprintf('\nNote that %s uses the parallel toolbox for acceleration if available, but otherwise may be slow!',mpm_params.R2s_fit_method)];
        case {'nlls_ols','nlls_wls1','nlls_wls2','nlls_wls3'}
            outstring=[outstring sprintf('\nWARNING that %s is very slow even if the parallel toolbox is available for acceleration! Consider using "wls1" instead (Edwards et al. ISMRM 2022).',mpm_params.R2s_fit_method)];
        otherwise
            outstring=[outstring sprintf('\nWARNING: But %s could not be found!',mpm_params.R2s_fit_method)];  
    end
    hmri_log(outstring,mpm_params.defflags);
end

% check whether there are enough echoes (neco4R2sfit) to estimate R2*
% (1) for basic R2* estimation, check only PDw images
mpm_params.basicR2s = false;
if mpm_params.PDwidx
    if length(mpm_params.input(mpm_params.PDwidx).TE) > (mpm_params.neco4R2sfit-1)
        mpm_params.basicR2s = true;
    end
end
% (2) for ESTATICS R2* estimation, check which (if any) contrast is usable
mpm_params.estaticsR2s = false*ones(1,mpm_params.ncon);
if mpm_params.proc.R2sOLS
    for ccon = 1:mpm_params.ncon
        mpm_params.estaticsR2s(ccon) = length(mpm_params.input(ccon).TE)>(mpm_params.neco4R2sfit-1);
    end
end
    
% if T2scorr enabled, must check there will be a R2s map generated!
if mpm_params.proc.PD.T2scorr && ~(any(mpm_params.estaticsR2s)||mpm_params.basicR2s)
    hmri_log(sprintf(['WARNING: not enough echoes available (minimum is %d) to ' ...
        '\ncalculate R2* map. No T2*-weighting bias correction can be ' ...
        '\napplied (T2scorr option). T2scorr disabled.']),mpm_params.defflags);
    mpm_params.proc.PD.T2scorr = 0;
end

% consistency check for number of echoes averaged for A calculation:
if mpm_params.PDwidx && mpm_params.T1widx 
    if mpm_params.proc.PD.nr_echoes_forA > size(mpm_params.input(mpm_params.T1widx).fnam,1)
        hmri_log(sprintf(['WARNING: number of T1w echoes to be averaged for PD calculation (%d)' ...
            '\nis bigger than the available number of echoes (%d). Setting nr_echoes_forA' ...
            '\nto %d, the maximum number of echoes available.'], mpm_params.proc.PD.nr_echoes_forA, ...
            size(mpm_params.input(mpm_params.T1widx).fnam,1), ...
            size(mpm_params.input(mpm_params.T1widx).fnam,1)),mpm_params.defflags);
        mpm_params.proc.PD.nr_echoes_forA = size(mpm_params.input(mpm_params.T1widx).fnam,1);
    end
end

% coregistration of all images to the PDw average (or TE=0 fit):
mpm_params.coreg = hmri_get_defaults('coreg2PDw');

% coregistration flags for weighted images
mpm_params.coreg_flags = hmri_get_defaults('coreg_flags');
hmri_log(sprintf('=== Registration Settings for weighted images ==='),mpm_params.nopuflags);
hmri_log(sprintf('Method: %s, Sampling: %s, Smoothing: %s', mpm_params.coreg_flags.cost_fun, mat2str(mpm_params.coreg_flags.sep), mat2str(mpm_params.coreg_flags.fwhm)),mpm_params.nopuflags);

% coregistration flags for B1 to PDw
mpm_params.coreg_bias_flags = hmri_get_defaults('coreg_bias_flags');
hmri_log(sprintf('=== Registration Settings for B1 bias images ==='),mpm_params.nopuflags);
hmri_log(sprintf('Method: %s, Sampling: %s, Smoothing: %s', mpm_params.coreg_bias_flags.cost_fun, mat2str(mpm_params.coreg_bias_flags.sep), mat2str(mpm_params.coreg_bias_flags.fwhm)),mpm_params.nopuflags);

% Prepare output for R1, PD, MT and R2* maps
RFsenscorr = fieldnames(mpm_params.proc.RFsenscorr);
B1transcorr = fieldnames(jobsubj.b1_type);
coutput = 0;

% R1 output: requires PDw and T1w input and optional B1 transmit (or
% UNICORT) and RF sensitivity maps for bias correction:
if (mpm_params.T1widx && mpm_params.PDwidx)
    coutput = coutput+1;
    mpm_params.qR1 = coutput;
    mpm_params.output(coutput).suffix = 'R1';
    mpm_params.output(coutput).units = 's-1';  
    mpm_params.output(coutput).descrip{1} = 'R1 map [s-1]';
    if mpm_params.proc.ISC.enabled
        mpm_params.output(coutput).descrip{end+1} = '- Imperfect Spoiling Correction applied';
    else
        mpm_params.output(coutput).descrip{end+1} = '- no Imperfect Spoiling Correction applied';
    end        
    switch B1transcorr{1}
        case 'no_B1_correction'
            mpm_params.output(coutput).descrip{end+1} = '- no B1+ bias correction applied';
        case 'UNICORT'
            mpm_params.output(coutput).descrip{end+1} = '- B1+ bias correction using UNICORT';
        otherwise
            mpm_params.output(coutput).descrip{end+1} = sprintf('- B1+ bias correction using provided B1 map (%s)',B1transcorr{1});
    end
    switch RFsenscorr{1}
        case 'RF_none'
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_us' % Unified Segmentation only applies to PD maps
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_once'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a single sensitivity measurement';
        case 'RF_per_contrast'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a per-contrast sensitivity measurement';
    end
else
    mpm_params.qR1 = 0;
end

% PD output: requires PDw and T1w input:
if (mpm_params.T1widx && mpm_params.PDwidx)
    coutput = coutput+1;
    mpm_params.qPD = coutput;
    if mpm_params.proc.PD.calibr && ...
            ~isfield(mpm_params.proc.RFsenscorr,'RF_none') && ...
            (~isempty(jobsubj.b1_trans_input)|| mpm_params.UNICORT.PD)
        mpm_params.output(coutput).suffix = 'PD';
        mpm_params.output(coutput).descrip{1} = 'PD map (water concentration) [p.u.]';
        mpm_params.output(coutput).descrip{end+1} = sprintf('- WM calibration (%g%%)', mpm_params.proc.PD.WMval);
        mpm_params.output(coutput).units = 'p.u.';
    else
        mpm_params.output(coutput).suffix = 'A';
        mpm_params.output(coutput).descrip{1} = 'A map (signal amplitude) [a.u.]';
        mpm_params.output(coutput).descrip{end+1} = '- no WM calibration';
        mpm_params.output(coutput).units = 'a.u.';
    end
    switch B1transcorr{1}
        case 'no_B1_correction'
            mpm_params.output(coutput).descrip{end+1} = '- no B1+ bias correction applied';
        case 'UNICORT'
            if mpm_params.UNICORT.PD
                mpm_params.output(coutput).descrip{end+1} = '- B1+ bias correction: UNICORT-estimated R1 map + UNICORT(R1)-estimated B1 map';
            else
                mpm_params.output(coutput).descrip{end+1} = '- B1+ bias correction: UNICORT-estimated R1 map used for A map calculation';
            end
        otherwise
            mpm_params.output(coutput).descrip{end+1} = sprintf('- B1+ bias correction using provided B1 map (%s)',B1transcorr{1});
    end
    switch RFsenscorr{1}
        case 'RF_none'
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_us' % Unified Segmentation only applies to PD maps
            if mpm_params.proc.PD.calibr
                mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on Unified Segmentation';
            else
                mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
            end
        case 'RF_once'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a single sensitivity measurement';
        case 'RF_per_contrast'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a per-contrast sensitivity measurement';
    end 
    if mpm_params.fullOLS 
        mpm_params.output(coutput).descrip{end+1} = '- R2* bias accounted for by TE=0 extrapolation (fullOLS option)';
    elseif mpm_params.proc.PD.T2scorr
        mpm_params.output(coutput).descrip{end+1} = '- R2* bias corrected for (T2scorr option)';            
    else
        mpm_params.output(coutput).descrip{end+1} = '- R2* bias not corrected for';            
    end
else
    mpm_params.qPD = 0;
end

% MT output: requires all three contrasts
if (mpm_params.T1widx && mpm_params.PDwidx && mpm_params.MTwidx)
    coutput = coutput+1;
    mpm_params.qMT = coutput;
    mpm_params.output(coutput).suffix = 'MTsat';
    mpm_params.output(coutput).descrip{1} = 'MT saturation map [p.u.]';
    mpm_params.output(coutput).units = 'p.u.';
    switch B1transcorr{1}
        case 'no_B1_correction'
            mpm_params.output(coutput).descrip{end+1} = '- no B1+ bias correction applied';
        case 'UNICORT'
            if mpm_params.UNICORT.MT
                mpm_params.output(coutput).descrip{end+1} = '- B1+ bias correction: UNICORT(R1)-estimated B1 map';
            else
                mpm_params.output(coutput).descrip{end+1} = '- no B1+ bias correction applied';
            end
        otherwise
           mpm_params.output(coutput).descrip{end+1} = sprintf('- B1+ bias correction using provided B1 map (%s)',B1transcorr{1});
    end
    switch RFsenscorr{1}
        case 'RF_none'
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_us' % Unified Segmentation only applies to PD maps
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_once'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a single sensitivity measurement';
        case 'RF_per_contrast'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a per-contrast sensitivity measurement';
    end

    % B1-correction method
    mpm_params.MTsatB1CorrectionModel = hmri_get_defaults('MTsatB1CorrectionModel');
    switch mpm_params.MTsatB1CorrectionModel
        case 'helms'
            mpm_params.MTsatB1CorrectionC = hmri_get_defaults('MTsatB1CorrectionHelmsC');
        case 'lipp'
            mpm_params.MTsatB1CorrectionC = hmri_get_defaults('MTsatB1CorrectionLippC');
        otherwise
            error('unknown MTsat B1 correction model ''%s''. Allowed models are ''helms'' and ''lipp''', mpm_params.MTsatB1CorrectionModel)
    end
    hmri_log(sprintf('INFO: Using MTsat B1 correction model ''%s'' with C = %g.', mpm_params.MTsatB1CorrectionModel, mpm_params.MTsatB1CorrectionC), mpm_params.defflags);
else
    mpm_params.qMT = 0;
end

if (mpm_params.PDwidx && mpm_params.MTwidx)
    if (mpm_params.input(mpm_params.MTwidx).TR == mpm_params.input(mpm_params.PDwidx).TR) ...
        && (mpm_params.input(mpm_params.MTwidx).fa == mpm_params.input(mpm_params.PDwidx).fa) % additional MTR image...
        coutput = coutput+1;
        mpm_params.qMTR = coutput;
        mpm_params.output(coutput).suffix = 'MTR';
        mpm_params.output(coutput).descrip{1} = 'Classic MTR image [a.u.] - not currently saved in Results';
        mpm_params.output(coutput).units = 'a.u.';
    else
        mpm_params.qMTR = 0;
    end
else
    mpm_params.qMTR = 0;
end

% R2* map: requires a minimum of 4 echoes with PDw contrast
if size(mpm_params.input(mpm_params.PDwidx).fnam,1) > (mpm_params.neco4R2sfit-1) 
    coutput = coutput+1;
    mpm_params.qR2s = coutput;
    mpm_params.output(coutput).suffix = 'R2s';
    mpm_params.output(coutput).descrip{1} = 'R2* map [s-1]';
    mpm_params.output(coutput).units = '[s-1]';
    if mpm_params.proc.R2sOLS
        mpm_params.output(coutput).descrip{end+1} = '- ESTATICS model (R2* map calculation)';
    end   
    mpm_params.output(coutput).descrip{end+1} = '- B1+ bias correction does not apply';
    switch RFsenscorr{1}
        case 'RF_none'
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_us' % Unified Segmentation only applies to PD maps
            mpm_params.output(coutput).descrip{end+1} = '- no RF sensitivity bias correction';
        case 'RF_once'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a single sensitivity measurement';
        case 'RF_per_contrast'
            mpm_params.output(coutput).descrip{end+1} = '- RF sensitivity bias correction based on a per-contrast sensitivity measurement';
    end     
else
    mpm_params.qR2s = 0;
end

% Summary of the output generated:
LogMsg = 'SUMMARY OF THE MAPS CALCULATED';
for coutput = 1:length(mpm_params.output)
    LogMsg = sprintf('%s\n\n(%d) %s (%s)',LogMsg, coutput, mpm_params.output(coutput).descrip{1} , mpm_params.output(coutput).suffix);
    for cdescrip = 2:length(mpm_params.output(coutput).descrip)
        LogMsg = sprintf('%s\n %s',LogMsg, mpm_params.output(coutput).descrip{cdescrip});
    end
end
hmri_log(LogMsg, mpm_params.nopuflags);

end


%% =======================================================================%
% To arrange the metadata structure for MPM calculation output.
%=========================================================================%
function metastruc = init_mpm_output_metadata(input_files, mpm_params)

proc.descrip = ['hMRI toolbox - ' mfilename '.m - map calculation'];
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
