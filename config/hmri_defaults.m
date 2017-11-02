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

% cleanup temporary directories. If set to true, all  
hmri_def.cleanup = true;
% settings for JSON metadata
hmri_def.json = struct('extended',false,'separate',true,'anonym','none',...
    'overwrite',true, 'indent','\t'); 
% recommended TPM for segmentation and spatial processing
hmri_def.TPM = fullfile(fileparts(fileparts(mfilename('fullpath'))),'etpm','eTPM.nii');
% default template for auto-reorientation
hmri_def.autoreorient_template = {fullfile(spm('dir'),'canonical','avg152T1.nii')};

%==========================================================================
% Default parameters for segmentation
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
% Ordinary Least Squares & fit at TE=0
%--------------------------------------------------------------------------
% create an Ordinary Least Squares R2* map?
hmri_def.R2sOLS = 1; 

% Define a coherent interpolation factor used all through the map creation
% process. Default is 3, but if you want to keep SNR and resolution as far
% as possible the same, it is recommended to use sinc interpolation (at
% least -4, in Siawoosh's experience -7 gives decent results)
hmri_def.interp = 3;

% Define the OLS fit as default. OLS fit at TE=0 is used instead of
% averaged contrast images for the map calculation
hmri_def.fullOLS = false;

%--------------------------------------------------------------------------
% PD maps processing parameters
%--------------------------------------------------------------------------
hmri_def.PDproc.PDmap    = 1;    % Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
hmri_def.PDproc.WBMaskTh = 0.1;  % Threshold for calculation of whole-brain mask from TPMs
hmri_def.PDproc.WMMaskTh = 0.95; % Threshold for calculation of white-matter mask from TPMs
hmri_def.PDproc.biasreg  = 10^(-5);
hmri_def.PDproc.biasfwhm = 50;
hmri_def.PDproc.nr_echoes_forA = 1; % NOTE: in order to minimize R2* bias 
    % on the PD estimates and gain in robustness for bias-field
    % correction, the first echo of the T1w series is used ("average"
    % calculated over the first echo only) for PD calculation. An average
    % over more echoes might be preferable when PD map's SNR is too poor,
    % but be aware that the gain in SNR will be balanced by an increased
    % R2* bias in PD values (in particular in the GM).
                         
%--------------------------------------------------------------------------
% RF sensitivity bias correction
%--------------------------------------------------------------------------
hmri_def.RFsens.smooth_kernel = 12;

%--------------------------------------------------------------------------
% quantitative maps: quality evaluation and realignment to MNI
%--------------------------------------------------------------------------
% creates a matlab structure containing markers of data quality
hmri_def.qMRI_maps.QA          = 1; 
% realigns qMRI maps to MNI: the following parameter corresponds to the
% realignment implemented as part of the map calculation (see
% hmri_create_MTProt.m). Left here for backward compatibility while it is
% recommended to rather reorient all images prior any processing using the
% Auto-Reorient module provided with the toolbox (type "help
% hmri_autoreorient" for details or open the SPM > Tools > hMRI Tools >
% Auto-Reorient module in the Batch GUI).
hmri_def.qMRI_maps.ACPCrealign = 0; 

%--------------------------------------------------------------------------
% Threshold values for qMRI maps
%--------------------------------------------------------------------------
hmri_def.qMRI_maps_thresh.R1       = 2000;
hmri_def.qMRI_maps_thresh.A        = 10^5;
hmri_def.qMRI_maps_thresh.R2s      = 10;
hmri_def.qMRI_maps_thresh.MTR      = 50;
hmri_def.qMRI_maps_thresh.MTR_synt = 50;
hmri_def.qMRI_maps_thresh.MT       = 5; 

%--------------------------------------------------------------------------
% MPM acquisition parameters and RF spoiling correction parameters
%--------------------------------------------------------------------------
% these value are initialised with defaults (v2k protocol - Prisma) for the
% first pass through this script only. They're updated at run-time with
% actual acquisition values (see hmri_MTProt.m).
% The coefficients for R.Deichmann steady state correction are also
% determined and stored as part of the defaults parameters for the current
% processing so they can be stored for reccord.

% Default MPMacq values to begin with
hmri_def.MPMacq.TE_mtw = 2.34;
hmri_def.MPMacq.TE_t1w = 2.34;
hmri_def.MPMacq.TE_pdw = 2.34;
hmri_def.MPMacq.TR_mtw = 24.5;
hmri_def.MPMacq.TR_t1w = 24.5;   % <-
hmri_def.MPMacq.TR_pdw = 24.5;   % <-
hmri_def.MPMacq.fa_mtw = 6;
hmri_def.MPMacq.fa_t1w = 21;     % <-
hmri_def.MPMacq.fa_pdw = 6;      % <-
hmri_def.MPMacq.tag    = 'v2k';

% Defining the MPMacq paramters distinguishing the different protocols
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

% Defining the RFCorr parameters for the different protocols
%--------------------------------------------------------------------------
% Antoine Lutti 15/01/09
% Correction parameters used in hmri_MTProt to correct for imperfect RF
% spoiling when a B1 map is loaded. Correction based on Preibisch and
% Deichmann's paper MRM 61:125-135 (2009). The values for P2_a and P2_b
% below were obtained using the code supplied by R. Deichmann with the
% experimental parameters used to get our PDw and T1w images. Correction
% parameters were calculated for the following parameter sets using
% T2 = 64 ms at 3T.
%
% 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
hmri_def.rfcorr.ClassicFIL.tag = 'Classic FIL protocol';
hmri_def.rfcorr.ClassicFIL.P2_a = [78.9228195006542,-101.113338489192,47.8783287525126];
hmri_def.rfcorr.ClassicFIL.P2_b = [-0.147476233142129,0.126487385091045,0.956824374979504];
hmri_def.rfcorr.ClassicFIL.RFCorr = true;
% 2) new FIL/Helms protocol
hmri_def.rfcorr.NewFILHelms.tag = 'New FIL/Helms protocol';
hmri_def.rfcorr.NewFILHelms.P2_a = [93.455034845930480,-120.5752858196904,55.911077913369060];
hmri_def.rfcorr.NewFILHelms.P2_b = [-0.167301931434861,0.113507432776106,0.961765216743606];
hmri_def.rfcorr.NewFILHelms.RFCorr = true;
% 3) Siemens product sequence protocol used in Lausanne (G Krueger)
hmri_def.rfcorr.SiemPrLausGK.tag = 'Siemens product Lausanne (GK) protocol';
hmri_def.rfcorr.SiemPrLausGK.P2_a = [67.023102027100880,-86.834117103841540,43.815818592349870];
hmri_def.rfcorr.SiemPrLausGK.P2_b = [-0.130876849571103,0.117721807209409,0.959180058389875];
hmri_def.rfcorr.SiemPrLausGK.RFCorr = true;
% 4) High-res (0.8mm) FIL protocol:
hmri_def.rfcorr.HResFIL.tag = 'High-res FIL protocol';
hmri_def.rfcorr.HResFIL.P2_a = [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
hmri_def.rfcorr.HResFIL.P2_b = [-0.218804328507184,0.178745853134922,0.939514554747592];
hmri_def.rfcorr.HResFIL.RFCorr = true;
% 5)NEW  High-res (0.8mm) FIL protocol:
hmri_def.rfcorr.NHResFIL.tag = 'New High-res FIL protocol';
hmri_def.rfcorr.NHResFIL.P2_a = [88.8623036106612,-114.526218941363,53.8168602253166];
hmri_def.rfcorr.NHResFIL.P2_b = [-0.132904017579521,0.113959390779008,0.960799295622202];
hmri_def.rfcorr.NHResFIL.RFCorr = true;
% 6)NEW  1mm protocol - seq version v2k:
hmri_def.rfcorr.v2k.tag = 'v2k protocol';
hmri_def.rfcorr.v2k.P2_a = [71.2817617982844,-92.2992876164017,45.8278193851731];
hmri_def.rfcorr.v2k.P2_b = [-0.137859046784839,0.122423212397157,0.957642744668469];
hmri_def.rfcorr.v2k.RFCorr = true;
% 7) 800um protocol - seq version v3* released used by MEG group:
hmri_def.rfcorr.v3star.tag = 'v3star protocol';
hmri_def.rfcorr.v3star.P2_a = [57.427573706259864,-79.300742898810441,39.218584751863879];
hmri_def.rfcorr.v3star.P2_b = [-0.121114060111119,0.121684347499374,0.955987357483519];
hmri_def.rfcorr.v3star.RFCorr = true;
% Unknwon protocol
hmri_def.rfcorr.Unknown.tag = 'Unknown protocol. No spoiling correction defined.';
hmri_def.rfcorr.Unknown.RFCorr = false;

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

end
