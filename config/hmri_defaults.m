function hmri_defaults
% Sets the defaults parameters which are used by the hMRI toolbox.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% DON'T MODIFY THIS FILE, IT CONTAINS THE REFERENCE DEFAULTS PARAMETERS.
% Please refer to hMRI-Toolbox\config\local\hmri_local_defaults.m to
% customise defaults parameters.  
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% FORMAT hmri_defaults
%__________________________________________________________________________
%
% THIS FILE SHOULD NOT BE MODIFIED. 
% To customize the hMRI-Toolbox defaults parameters so they match your own
% site- or protocol-specific setup, please refer to the defaults files in
% hMRI-Toolbox\config\local. In particular, use "hmri_local_defaults.m".
% Make a copy with meaningful name, modify as desired and select as general
% defaults file in the "Configure toolbox" branch of the hMRI-Toolbox.
%
% The structure and content of this file are largely inspired by the
% equivalent file in SPM.
%
% DOCUMENTATION
% A brief description of each parameter is provided together with
% guidelines and recommendations to modify these parameters. With few
% exceptions, parameters should ONLY be MODIFIED and customized BY ADVANCED
% USERS, having a good knowledge of the underlying algorithms and
% implementation. Please refer to the documentation in the github WIKI for
% more details.  
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liege, Belgium

% Global hmri_def variable used across the whole toolbox
global hmri_def

% Specifying the research centre - to be customized in the local
% configuration file (config/local/hmri_local_defaults.m). Not mandatory.
hmri_def.centre = 'centre' ; 

% Defaults customised defaults file location
hmri_def.local_defaults = {fullfile(fileparts(mfilename('fullpath')),'local','hmri_local_defaults.m')};

%==========================================================================
% Common processing parameters 
%==========================================================================

% cleanup temporary directories. If set to true, all temporary directories
% are deleted at the end of map creation, only the "Results" directory and
% "Supplementary" subdirectory are kept. Setting "cleanup" to "false" might
% be convenient if one desires to have a closer look at intermediate
% processing steps. Otherwise "cleanup = true" is recommended for saving
% disk space.
hmri_def.cleanup = true;
% settings for JSON metadata: by default, separate JSON files are used to
% store the metadata (information on data acquisition and processing,
% tracking of input and output files), as JSON-formatted, tab-indented
% text. The following settings are recommended. No modification currently
% foreseen as useful...
hmri_def.json = struct('extended',false,'separate',true,'anonym','none',...
    'overwrite',true, 'indent','\t'); 
% recommended TPM for segmentation and spatial processing. The hMRI toolbox
% provides a series of tissue probability maps. These TPMs could be
% replaced by other TPMs, to better match the population studied. 
% ADVANCED USER ONLY.
hmri_def.TPM = fullfile(fileparts(fileparts(which('hmri_defaults.m'))),'etpm','eTPM.nii');
% default template for auto-reorientation. The template can be selected
% within the Auto-reorient module. The following is the default suggested
% for T1w images. Please refer to the Auto-reorient documentation for an
% appropriate choice of the template.
hmri_def.autoreorient_template = {fullfile(spm('dir'),'canonical','avg152T1.nii')};

%==========================================================================
% Default parameters for segmentation
% ADVANCED USERS ONLY!
% hmri_def.segment is effectively the job to be handed to spm_preproc_run
% By default, parameters are set to
% - create tissue class images (c*) in the native space of the source image
%   (tissue(*).native = [1 0]) for tissue classes 1-5
% - save both BiasField and BiasCorrected volume (channel.write = [1 1])
% - recommended values from SPM12 (October 2017)
%==========================================================================

% hmri_def.segment.channel.vols = cell array of file names, 
%                       must be defined before calling spm_preproc_run
hmri_def.segment.channel.biasreg = 0.001;
hmri_def.segment.channel.biasfwhm = 60;
% hmri_def.segment.channel.write = [0 0]; % save nothing
% hmri_def.segment.channel.write = [1 0]; % save BiasField
% hmri_def.segment.channel.write = [0 1]; % save BiasCorrected volume
hmri_def.segment.channel.write = [1 1]; % save BiasField and BiasCorrected volume

hmri_def.segment.tissue(1).tpm = {[hmri_def.TPM ',1']};
hmri_def.segment.tissue(1).ngaus = 1;
hmri_def.segment.tissue(1).native = [1 0];
hmri_def.segment.tissue(1).warped = [0 0];
hmri_def.segment.tissue(2).tpm = {[hmri_def.TPM ',2']};
hmri_def.segment.tissue(2).ngaus = 1;
hmri_def.segment.tissue(2).native = [1 0];
hmri_def.segment.tissue(2).warped = [0 0];
hmri_def.segment.tissue(3).tpm = {[hmri_def.TPM ',3']};
hmri_def.segment.tissue(3).ngaus = 2;
hmri_def.segment.tissue(3).native = [1 0];
hmri_def.segment.tissue(3).warped = [0 0];
hmri_def.segment.tissue(4).tpm = {[hmri_def.TPM ',4']};
hmri_def.segment.tissue(4).ngaus = 3;
hmri_def.segment.tissue(4).native = [1 0];
hmri_def.segment.tissue(4).warped = [0 0];
hmri_def.segment.tissue(5).tpm = {[hmri_def.TPM ',5']};
hmri_def.segment.tissue(5).ngaus = 4;
hmri_def.segment.tissue(5).native = [1 0];
hmri_def.segment.tissue(5).warped = [0 0];
hmri_def.segment.tissue(6).tpm = {[hmri_def.TPM ',6']};
hmri_def.segment.tissue(6).ngaus = 2;
hmri_def.segment.tissue(6).native = [0 0];
hmri_def.segment.tissue(6).warped = [0 0];
hmri_def.segment.warp.mrf = 1;
hmri_def.segment.warp.cleanup = 1;
hmri_def.segment.warp.reg = [0 0.001 0.5 0.05 0.2];
hmri_def.segment.warp.affreg = 'mni';
hmri_def.segment.warp.fwhm = 0;
hmri_def.segment.warp.samp = 3;
hmri_def.segment.warp.write = [0 0];

%==========================================================================
% R1/PD/R2s/MT map creation parameters
%==========================================================================

%--------------------------------------------------------------------------
% Coregistration of all input images to the average (or TE=0 fit) PDw image
%--------------------------------------------------------------------------
% The coregistration step can be disabled using the following flag (not
% recommended). ADVANCED USER ONLY. 
hmri_def.coreg2PDw = 1; 

%--------------------------------------------------------------------------
% Coregistration parameters for weighted images to the average 
% (or TE=0 fit) PDw image; see spm_coreg.m for details
%--------------------------------------------------------------------------
% For high resolution (< 0.8 mm) data, we have found that adding extra 
% steps to the registration can improve map quality, e.g. for 0.5 mm data
% we have found
%   hmri_def.coreg_flags.sep=[4 2 1 0.5];
% to be needed to avoid registration artefacts, but this comes at the 
% expense of slowing down map creation
hmri_def.coreg_flags.sep=[4 2];
hmri_def.coreg_flags.cost_fun = 'nmi';
hmri_def.coreg_flags.fwhm = [7 7];

%--------------------------------------------------------------------------
% Coregistration parameters for B1 bias maps to the average (or TE=0 fit) 
% PDw image; see spm_coreg.m for details
%--------------------------------------------------------------------------
hmri_def.coreg_bias_flags.sep=[4 2];
hmri_def.coreg_bias_flags.cost_fun = 'nmi';
hmri_def.coreg_bias_flags.fwhm = [7 7];

%--------------------------------------------------------------------------
% Decide whether to use small angle approximation when computing R1 and PD
%--------------------------------------------------------------------------
hmri_def.small_angle_approx = true;

%--------------------------------------------------------------------------
% Ordinary Least Squares & fit at TE=0
%--------------------------------------------------------------------------
% Create a Least Squares R2* map. The ESTATICS model is applied
% to calculate a common R2* map from all available contrasts. 
% ADVANCED USER ONLY.
hmri_def.R2sOLS = 1; 

% Choose method of R2* fitting; this is a trade-off between speed and
% accuracy.
% - 'OLS' (classic ESTATICS ordinary least squares model; fast but not as 
%         accurate as WLS1)
% - 'WLS1' (weighted least squares with one iteration; slower than OLS but 
%          significantly more accurate. Uses parallelization over voxels to 
%          speed up calculation)
% - 'WLS2', 'WLS3' (weighted least squares with 2 or 3 iterations; did not
%                   show great improvement over WLS1 in test data)
% - 'NLLS_OLS', 'NLLS_WLS1' (nonlinear least squares using either OLS or
%                            WLS1 to provide the initial parameter
%                            estimates; very slow, even with
%                            parallelization over voxels. Recommend WLS1
%                            instead)
hmri_def.R2s_fit_method = 'OLS';

% Minimum number of echoes to calculate R2s map. Strictly speaking, the
% minimum is 2. For a robust estimation, the minimum number of echoes
% required depends on many factors, amongst which: 
% - SNR/resolution
% - distribution/spacing between TEs: note that early echoes are more
%   affected by the specific contrast, violating the assumption of a common
%   decay between contrasts. 
% - number of contrasts available (fewer echoes per contrast required for 3
%   (PDw, T1w, MTw) contrasts as compared to 2 or even 1)
% To be on the safe side, a minimum of 6 echoes is recommended (ESTATICS
% paper). Further studies are required to come up with more detailed and
% informed guidelines. Use fewer echoes at your own risk...! 
hmri_def.neco4R2sfit = 4;

% Define a coherent interpolation factor used all through the map creation
% process. Default is 3, but if you want to keep SNR and resolution as far
% as possible the same, it is recommended to use sinc interpolation (at
% least -4, in Siawoosh's experience -7 gives decent results). 
% ADVANCED USER ONLY. 
hmri_def.interp = 3;

% Define the OLS fit as default. OLS fit at TE=0 for each contrast is used
% instead of averaged contrast images for the map calculation.
% ADVANCED USER ONLY. 
hmri_def.fullOLS = true;

%--------------------------------------------------------------------------
% Usage of UNICORT-derived B1 maps for PD and/or MT maps calculation
% ADVANCED USER ONLY.
% WARNING: this method has not been validated for PD and MT calculation!
%--------------------------------------------------------------------------
hmri_def.UNICORT.PD = false;
hmri_def.UNICORT.MT = false;

%--------------------------------------------------------------------------
% PD maps processing parameters
% ADVANCED USER ONLY.
%--------------------------------------------------------------------------
hmri_def.PDproc.calibr    = 1;   % Calibration of the PD map (if PDw, T1w, 
    % B1map available and RF sensitivity bias correction applied somehow)
    % based on PD(WM) = 69% [Tofts 2003]. 
hmri_def.PDproc.WBMaskTh = 0.1;  % Threshold for calculation of whole-brain mask from TPMs
hmri_def.PDproc.WMMaskTh = 0.95; % Threshold for calculation of white-matter mask from TPMs
hmri_def.PDproc.biasreg  = 10^(-5);
hmri_def.PDproc.biasfwhm = 50;
hmri_def.PDproc.nr_echoes_forA = 6; % NOTE: in order to minimize R2* bias 
    % on the PD estimates and gain in robustness for bias-field
    % correction, the number of echoes should be minimum ("average"
    % calculated over the first echo only) for PD calculation. However,
    % with T2*-weighting bias correction (see below), a higher number of
    % echoes is preferred in order to provide good SNR. Note that when
    % "fullOLS = true", this parameter has no impact whatsovever.
hmri_def.PDproc.T2scorr = 1; % to correct A map for T2*-weighting bias 
    % before PD map calculation. The correction is not required when
    % "fullOLS = true" (TE=0 fit) and will be automatically disabled.

%--------------------------------------------------------------------------
% RF sensitivity bias correction: kernel for sensitivity map smoothing.
% ADVANCED USER ONLY.
%--------------------------------------------------------------------------
hmri_def.RFsens.smooth_kernel = 12;

%--------------------------------------------------------------------------
% quantitative maps: quality evaluation and realignment to MNI
%--------------------------------------------------------------------------
% creates a matlab structure containing markers of data quality
hmri_def.qMRI_maps.QA          = 1; 
% realigns qMRI maps to MNI: the following parameter corresponds to the
% realignment implemented as part of the map calculation (see
% hmri_create_MTProt.m). Left here for backward compatibility. It is
% STRONGLY RECOMMENDED to reorient all images prior any processing using 
% the Auto-Reorient module provided with the toolbox (type "help
% hmri_autoreorient" for details or open the SPM > Tools > hMRI Tools >
% Auto-Reorient module in the Batch GUI). ADVANCED USER ONLY.
hmri_def.qMRI_maps.ACPCrealign = 0; 

%--------------------------------------------------------------------------
% Threshold values for qMRI maps
% The thresholds are meant to discard outliers generally due to low SNR in
% some brain areas, leading to physical non-sense values. Thresholding is
% required to process further the maps generated, when e.g. used
% segmentation algorithms make assumptions incompatible with existing
% outliers.
% NOTE that thresholding modifies the signal distribution and may alter
% the statistical results.
% ADVANCED USER ONLY.
%--------------------------------------------------------------------------
hmri_def.qMRI_maps_thresh.R1       = 2000; % 1000*[s-1]
hmri_def.qMRI_maps_thresh.A        = 10^5; % [a.u.] based on input images with intensities ranging approx. [0 4096].
hmri_def.qMRI_maps_thresh.R2s      = 10;   % 1000*[s-1]
hmri_def.qMRI_maps_thresh.MTR      = 50;
hmri_def.qMRI_maps_thresh.MTR_synt = 50;
hmri_def.qMRI_maps_thresh.MT       = 5;    % [p.u.]

%--------------------------------------------------------------------------
% MPM acquisition parameters and RF spoiling correction parameters
%--------------------------------------------------------------------------
% ACQUISITION PARAMETERS: these values are initialised with defaults (v2k
% protocol - Prisma) and are updated at run-time with actual acquisition
% values (see hmri_create_MTProt.m). If TR/TE/FA cannot be determined from
% the input images, the following values will be used. If they don't match
% your own protocol values and if no TR/TE/FA values can be retrieved by
% the toolbox from your data, the following values should be adapted in the
% local defaults file. 
% ADVANCED USER ONLY
hmri_def.MPMacq.TE_mtw = [2.34:2.34:14.04]'; %#ok<*NBRAK> % [ms] defined as column vector!
hmri_def.MPMacq.TE_t1w = [2.34:2.34:18.72]'; % [ms]
hmri_def.MPMacq.TE_pdw = [2.34:2.34:18.72]'; % [ms]
hmri_def.MPMacq.TR_mtw = 24.5; % [ms]
hmri_def.MPMacq.TR_t1w = 24.5; % [ms]
hmri_def.MPMacq.TR_pdw = 24.5; % [ms]
hmri_def.MPMacq.fa_mtw = 6;    % [deg]
hmri_def.MPMacq.fa_t1w = 21;   % [deg]
hmri_def.MPMacq.fa_pdw = 6;    % [deg]
hmri_def.MPMacq.tag    = 'v2k';

% IMPREFECT RF SPOILING CORRECTION PARAMETERS 
% (Preibisch and Deichmann, MRM 61:125-135 (2009))
% No run-time calculation of the correction parametes is currently
% performed in the toolbox. They've been calculated and stored here for a
% series of standard protocols, and are selected at run-time according to
% effective TR and FA values (more precisely: [TR_pdw TR_t1w fa_pdw
% fa_t1w]). If the used TR and FA values don't match any of the predefined
% types, no correction is applied. Additional protocol types and correction
% factors can be computed off-line (see more details below) and added to
% the local defaults file to apply proper RF correction.
% ADVANCED USER ONLY.

% A) Defining the MPMacq parameters distinguishing the different protocols
%--------------------------------------------------------------------------
% Using the following parameter order: [TR_pdw TR_t1w fa_pdw fa_t1w]
% NOTE: all tags MUST 
% - start with a letter, and 
% - include only letters, numbers or underscore, i.e. NO space.
% as these names are used to define a structure fieldname with the protocol 
% parameters.
%
% 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
% PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
hmri_def.MPMacq_set.names{1} = 'Classic FIL protocol';
hmri_def.MPMacq_set.tags{1}  = 'ClassicFIL';
hmri_def.MPMacq_set.vals{1}  = [23.7 18.7 6 20];
% 2) new FIL/Helms protocol
% PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
hmri_def.MPMacq_set.names{2} = 'New FIL/Helms protocol';
hmri_def.MPMacq_set.tags{2}  = 'NewFILHelms';
hmri_def.MPMacq_set.vals{2}  = [24.5 24.5 5 29];
% 3) Siemens product sequence protocol used in Lausanne (G Krueger)
% PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
hmri_def.MPMacq_set.names{3} = 'Siemens product Lausanne (GK) protocol';
hmri_def.MPMacq_set.tags{3}  = 'SiemPrLausGK';
hmri_def.MPMacq_set.vals{3}  = [24.0 19.0 6 20];
% 4) High-res (0.8mm) FIL protocol:
% PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
hmri_def.MPMacq_set.names{4} = 'High-res FIL protocol';
hmri_def.MPMacq_set.tags{4}  = 'HResFIL';
hmri_def.MPMacq_set.vals{4}  = [23.7 23.7 6 28];
% 5)NEW  High-res (0.8mm) FIL protocol:
% PD-weighted: TR=25.25ms; a=5deg; T1-weighted: TR=TR=25.25ms; a=29deg
hmri_def.MPMacq_set.names{5} = 'New High-res FIL protocol';
hmri_def.MPMacq_set.tags{5}  = 'NHResFIL';
hmri_def.MPMacq_set.vals{5}  = [25.25 25.25 5 29];
% 6)NEW  1mm protocol - seq version v2k:
% PD-weighted: TR=24.5ms; a=6deg; T1-weighted: TR=24.5ms; a=21deg
hmri_def.MPMacq_set.names{6} = 'v2k protocol';
hmri_def.MPMacq_set.tags{6}  = 'v2k';
hmri_def.MPMacq_set.vals{6}  = [24.5 24.5 6 21];
% 7) 800um protocol - seq version v3* released used by MEG group:
% TR = 25ms for all volumes; flipAngles = [6, 21 deg] for PDw and T1w
% Correction parameters below were determined via Bloch-Torrey 
% simulations but end result agrees well with EPG-derived correction 
% for this RF spoiling increment of 137 degrees.
% See: Callaghan et al. ISMRM, 2015, #1694
hmri_def.MPMacq_set.names{7} = 'v3star protocol';
hmri_def.MPMacq_set.tags{7}  = 'v3star';
hmri_def.MPMacq_set.vals{7}  = [25 25 6 21];

% B) Defining the imperfect spoiling correction parameters for the
% different protocols 
%--------------------------------------------------------------------------
% Antoine Lutti 15/01/09
% The following correction parameters are based on 
% [Preibisch and Deichmann, Magn Reson Med 61:125-135 (2009)]. The values
% for P2_a and P2_b below were obtained using the code supplied by R.
% Deichmann with the experimental parameters used for the standard MPM
% protocol and assuming T2 = 64 ms at 3T.

% Given the possible confusion and resulting mistake (imperfect spoiling
% correction applied to the wrong sequence), when TR and FA values match
% one of the listed cases below, the option is disabled by default. 
% When enabling the imperfect spoiling correction, make sure the
% coefficients retrieved in the list below are definitely calculated for
% the protocol used!
hmri_def.imperfectSpoilCorr.enabled = false;

% 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
hmri_def.imperfectSpoilCorr.ClassicFIL.tag = 'Classic FIL protocol';
hmri_def.imperfectSpoilCorr.ClassicFIL.P2_a = [78.9228195006542,-101.113338489192,47.8783287525126];
hmri_def.imperfectSpoilCorr.ClassicFIL.P2_b = [-0.147476233142129,0.126487385091045,0.956824374979504];
hmri_def.imperfectSpoilCorr.ClassicFIL.small_angle_approx = true;
% 2) new FIL/Helms protocol
hmri_def.imperfectSpoilCorr.NewFILHelms.tag = 'New FIL/Helms protocol';
hmri_def.imperfectSpoilCorr.NewFILHelms.P2_a = [93.455034845930480,-120.5752858196904,55.911077913369060];
hmri_def.imperfectSpoilCorr.NewFILHelms.P2_b = [-0.167301931434861,0.113507432776106,0.961765216743606];
hmri_def.imperfectSpoilCorr.NewFILHelms.small_angle_approx = true;
% 3) Siemens product sequence protocol used in Lausanne (G Krueger)
hmri_def.imperfectSpoilCorr.SiemPrLausGK.tag = 'Siemens product Lausanne (GK) protocol';
hmri_def.imperfectSpoilCorr.SiemPrLausGK.P2_a = [67.023102027100880,-86.834117103841540,43.815818592349870];
hmri_def.imperfectSpoilCorr.SiemPrLausGK.P2_b = [-0.130876849571103,0.117721807209409,0.959180058389875];
hmri_def.imperfectSpoilCorr.SiemPrLausGK.small_angle_approx = true;
% 4) High-res (0.8mm) FIL protocol:
hmri_def.imperfectSpoilCorr.HResFIL.tag = 'High-res FIL protocol';
hmri_def.imperfectSpoilCorr.HResFIL.P2_a = [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
hmri_def.imperfectSpoilCorr.HResFIL.P2_b = [-0.218804328507184,0.178745853134922,0.939514554747592];
hmri_def.imperfectSpoilCorr.HResFIL.small_angle_approx = true;
% 5)NEW  High-res (0.8mm) FIL protocol:
hmri_def.imperfectSpoilCorr.NHResFIL.tag = 'New High-res FIL protocol';
hmri_def.imperfectSpoilCorr.NHResFIL.P2_a = [88.8623036106612,-114.526218941363,53.8168602253166];
hmri_def.imperfectSpoilCorr.NHResFIL.P2_b = [-0.132904017579521,0.113959390779008,0.960799295622202];
hmri_def.imperfectSpoilCorr.NHResFIL.small_angle_approx = true;
% 6)NEW  1mm protocol - seq version v2k:
hmri_def.imperfectSpoilCorr.v2k.tag = 'v2k protocol';
hmri_def.imperfectSpoilCorr.v2k.P2_a = [71.2817617982844,-92.2992876164017,45.8278193851731];
hmri_def.imperfectSpoilCorr.v2k.P2_b = [-0.137859046784839,0.122423212397157,0.957642744668469];
hmri_def.imperfectSpoilCorr.v2k.small_angle_approx = true;
% 7) 800um protocol - seq version v3* released used by MEG group:
hmri_def.imperfectSpoilCorr.v3star.tag = 'v3star protocol';
hmri_def.imperfectSpoilCorr.v3star.P2_a = [57.427573706259864,-79.300742898810441,39.218584751863879];
hmri_def.imperfectSpoilCorr.v3star.P2_b = [-0.121114060111119,0.121684347499374,0.955987357483519];
hmri_def.imperfectSpoilCorr.v3star.small_angle_approx = true;
% Unknown protocol
hmri_def.imperfectSpoilCorr.Unknown.tag = 'Unknown protocol. No spoiling correction defined.';

%--------------------------------------------------------------------------
% B1 mapping processing parameters 
%--------------------------------------------------------------------------
% All defaults for B1 map calculation are defined in 
% hMRI-Toolbox\config\hmri_b1_standard_defaults.m. 
% See examples of local customization in the hMRI-Toolbox\config\local
% directory. 

hmri_b1_standard_defaults;
%==========================================================================
% Maps processing parameters
%==========================================================================

%--------------------------------------------------------------------------
% US segmentation parameters
%--------------------------------------------------------------------------

% recommended TPM for segmentation
hmri_def.proc.TPM = hmri_def.TPM ;
% Use the same as for the maps creation but one could (want to) use another
% one at some point. 
% Map creation works with "standard" weighted-MR images to build the 
% parametric maps. In the end these parametric maps taken together for a 
% multichannel-segmention could show more details (for example subcortical 
% nuclei?) and would therefore require a specific TPM. This TPM is of 
% course still to be built at the moment...

% Flags to write out posterior tissue classes in native & warped space
% - GM/WM/CSF -> write warped, mod+unmod, and native, native+dartelImp.
% - others -> nothing
hmri_def.proc.w_native = [[1 1];[1 1];[1 1];[0 0];[0 0];[0 0]];
hmri_def.proc.w_warped = [[1 1];[1 1];[1 1];[0 0];[0 0];[0 0]];
% Number of Gaussians per tissue class
hmri_def.proc.nGauss = [2 2 2 3 4 2]; % originally in SPM [1 1 2 3 4 2]

%==========================================================================
% Generation of error maps
%==========================================================================
% errormaps
hmri_def.errormaps  = false;

% lower bounds on R1, PD and MTsat to avoid dividing by small numbers
% in parameter-wise SNR estimates of R1/dR1, PD/dPD and MTsat/dMTsat
hmri_def.qMRI_maps_thresh.dR1 = 1e-4;
hmri_def.qMRI_maps_thresh.dPD = 1e-2;
hmri_def.qMRI_maps_thresh.dMT = 1e-4;

% upper bounds on SNR maps dR1/R1, dPD/PD and dMTsat/MTsat
hmri_def.qMRI_maps_thresh.SNR_R1= 1e3;
hmri_def.qMRI_maps_thresh.SNR_PD= 1e3;
hmri_def.qMRI_maps_thresh.SNR_MT= 1e3;

% weighted combination - see Mohammadi et al., NeuroImage, 2022 for details
% might have to be adjusted for different dataset
hmri_def.wcombparams.kt         = 10;   % This parameter is relevant and might have to be adjusted per protocol (defined in percent). 
% Examples can be found in Mohammadi et al., NeuroImage, 2022, Supplementary Materials: S1: Efficiency of robust combination and the Fermi function
hmri_def.wcombparams.am         = true; % True: arithmetic mean is written out, in addition to robust combination
hmri_def.wcombparams.errormaps  = true; % True: robust combination of error maps is generated and written out
% The following parameters are for experts only.
hmri_def.wcombparams.res        = -4;   % Resampling factor, determines interpolation method (details can be found in spm_slice_vol.m)
hmri_def.wcombparams.smthk      = 0;    % If larger than zero, the error maps will be spatially smoothed with the a Gaussian kernel of FWHM = smthk
hmri_def.wcombparams.dim        = 3;    % Defines dimention along which slices along which weightes will be determined

end
