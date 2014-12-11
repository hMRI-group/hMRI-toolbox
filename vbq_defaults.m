function vbq_defaults
% Sets the defaults which are used by the VBQ toolbox.
%
% FORMAT vbq_defaults
%_______________________________________________________________________
%
% This file can be customised to any the site/person own setup.
% Individual users can make copies which can be stored on their own
% matlab path. Make sure then that your 'vbq_defaults' is the first one
% found in the path. See matlab documentation for details on setting path.
%
% Care must be taken when modifying this file!
%
% The structure and content of this file are largely inspired by the
% equivalent file in SPM.
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: vbq_defaults.m 30 2013-11-27 14:50:20Z christophe $

%%
global vbq_def

%% %%%%%%%%%%%%%%%%%%%%% Global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specifying the lab
vbq_def.centre = 'crc' ; % can be 'fil' or 'lren' or 'crc'

%% %%%%%%%%%%%%%%%%% Common processing parameters %%%%%%%%%%%%%%%%%%%%%
% These parameters are either parameters that are fixed for all sites or
% recommended values. They can also be changed in a site-specific way at
% run-time.

%% Processing of PD maps
vbq_def.PDproc.PDmap    = 1;    % Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
vbq_def.PDproc.WBMaskTh = 0.1;  % Threshold for calculation of whole-brain mask from TPMs
vbq_def.PDproc.WMMaskTh = 0.95; % Threshold for calculation of white-matter mask from TPMs
vbq_def.PDproc.biasreg  = 10^(-5);
vbq_def.PDproc.biasfwhm = 50;

%% Threshold values for saving of the qMRI maps
vbq_def.qMRI_maps_thresh.R1       = 2000;
vbq_def.qMRI_maps_thresh.A        = 10^5;
vbq_def.qMRI_maps_thresh.R2s      = 10;
vbq_def.qMRI_maps_thresh.MTR      = 50;
vbq_def.qMRI_maps_thresh.MTR_synt = 50;
vbq_def.qMRI_maps_thresh.MT       = 5;

%% MPM acquisition parameters and RF spoiling correction parameters
% these value are initialised with defaults (v2k protocol - Prisma) for the
% first pass through this script only. They're updated at run-time with
% actual acquisition values (see vbq_MTProt.m).
% The coefficients for R.Deichmann steady state correction are also
% determined and stored as part of the defaults parameters for the current
% processing so they can be stored for reccord.

% Default MPMacq values to begin with
%---------------------------------------------------------------------
vbq_def.MPMacq.TE_mtw = 2.34;
vbq_def.MPMacq.TE_t1w = 2.34;
vbq_def.MPMacq.TE_pdw = 2.34;
vbq_def.MPMacq.TR_mtw = 24.5;
vbq_def.MPMacq.TR_t1w = 24.5;   % <-
vbq_def.MPMacq.TR_pdw = 24.5;   % <-
vbq_def.MPMacq.fa_mtw = 6;
vbq_def.MPMacq.fa_t1w = 21;     % <-
vbq_def.MPMacq.fa_pdw = 6;      % <-
vbq_def.MPMacq.tag    = 'v2k';

% Defining the MPMacq paramters distinguishing the different protocols
%---------------------------------------------------------------------
% Using the following parameter order: [TR_pdw TR_t1w fa_pdw fa_t1w]
% NOTE: 
% all the tags MUST 
% - start with a letter, and 
% - include only letters, numbers or underscore, i.e. NO space.
% as these names are used to define a structure fieldname with the protocol 
% parameters.
%
% 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
% PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
vbq_def.MPMacq_set.names{1} = 'Classic FIL protocol';
vbq_def.MPMacq_set.tags{1}  = 'ClassicFIL';
vbq_def.MPMacq_set.vals{1}  = [23.7 18.7 6 20];
% 2) new FIL/Helms protocol
% PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
vbq_def.MPMacq_set.names{2} = 'New FIL/Helms protocol';
vbq_def.MPMacq_set.tags{2}  = 'NewFILHelms';
vbq_def.MPMacq_set.vals{2}  = [24.5 24.5 5 29];
% 3) Siemens product sequence protocol used in Lausanne (G Krueger)
% PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
vbq_def.MPMacq_set.names{3} = 'Siemens product Lausanne (GK) protocol';
vbq_def.MPMacq_set.tags{3}  = 'SiemPrLausGK';
vbq_def.MPMacq_set.vals{3}  = [24.0 19.0 6 20];
% 4) High-res (0.8mm) FIL protocol:
% PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
vbq_def.MPMacq_set.names{4} = 'High-res FIL protocol';
vbq_def.MPMacq_set.tags{4}  = 'HResFIL';
vbq_def.MPMacq_set.vals{4}  = [23.7 23.7 6 28];
% 5)NEW  High-res (0.8mm) FIL protocol:
% PD-weighted: TR=25.25ms; a=5deg; T1-weighted: TR=TR=25.25ms; a=29deg
vbq_def.MPMacq_set.names{5} = 'New High-res FIL protocol';
vbq_def.MPMacq_set.tags{5}  = 'NHResFIL';
vbq_def.MPMacq_set.vals{5}  = [25.25 25.25 5 29];
% 6)NEW  1mm protocol - seq version v2k:
% PD-weighted: TR=24.5ms; a=6deg; T1-weighted: TR=24.5ms; a=21deg
vbq_def.MPMacq_set.names{6} = 'v2k protocol';
vbq_def.MPMacq_set.tags{6}  = 'v2k';
vbq_def.MPMacq_set.vals{6}  = [24.5 24.5 6 21];

% Defining the RFCorr parameters for the different protocols
%---------------------------------------------------------------------
% Antoine Lutti 15/01/09
% Correction parameters used in vbq_MTProt to correct for imperfect RF
% spoiling when a B1 map is loaded. Correction based on Preibisch and
% Deichmann's paper MRM 61:125-135 (2009). The values for P2_a and P2_b
% below were obtained using the code supplied by R. Deichmann with the
% experimentalparameters used to get our PDw and T1w images. Correction
% parameters were calculated for the following parameter sets using
% T2 = 64 ms at 3T.
%
% 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
vbq_def.rfcorr.ClassicFIL.tag = 'Classic FIL protocol';
vbq_def.rfcorr.ClassicFIL.P2_a = [78.9228195006542,-101.113338489192,47.8783287525126];
vbq_def.rfcorr.ClassicFIL.P2_b = [-0.147476233142129,0.126487385091045,0.956824374979504];
vbq_def.rfcorr.ClassicFIL.RFCorr = true;
% 2) new FIL/Helms protocol
vbq_def.rfcorr.NewFILHelms.tag = 'New FIL/Helms protocol';
vbq_def.rfcorr.NewFILHelms.P2_a = [93.455034845930480,-120.5752858196904,55.911077913369060];
vbq_def.rfcorr.NewFILHelms.P2_b = [-0.167301931434861,0.113507432776106,0.961765216743606];
vbq_def.rfcorr.NewFILHelms.RFCorr = true;
% 3) Siemens product sequence protocol used in Lausanne (G Krueger)
vbq_def.rfcorr.SiemPrLausGK.tag = 'Siemens product Lausanne (GK) protocol';
vbq_def.rfcorr.SiemPrLausGK.P2_a = [67.023102027100880,-86.834117103841540,43.815818592349870];
vbq_def.rfcorr.SiemPrLausGK.P2_b = [-0.130876849571103,0.117721807209409,0.959180058389875];
vbq_def.rfcorr.SiemPrLausGK.RFCorr = true;
% 4) High-res (0.8mm) FIL protocol:
vbq_def.rfcorr.HResFIL.tag = 'High-res FIL protocol';
vbq_def.rfcorr.HResFIL.P2_a = [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
vbq_def.rfcorr.HResFIL.P2_b = [-0.218804328507184,0.178745853134922,0.939514554747592];
vbq_def.rfcorr.HResFIL.RFCorr = true;
% 5)NEW  High-res (0.8mm) FIL protocol:
vbq_def.rfcorr.NHResFIL.tag = 'New High-res FIL protocol';
vbq_def.rfcorr.NHResFIL.P2_a = [88.8623036106612,-114.526218941363,53.8168602253166];
vbq_def.rfcorr.NHResFIL.P2_b = [-0.132904017579521,0.113959390779008,0.960799295622202];
vbq_def.rfcorr.NHResFIL.RFCorr = true;
% 6)NEW  1mm protocol - seq version v2k:
vbq_def.rfcorr.v2k.tag = 'v2k protocol';
vbq_def.rfcorr.v2k.P2_a = [71.2817617982844,-92.2992876164017,45.8278193851731];
vbq_def.rfcorr.v2k.P2_b = [-0.137859046784839,0.122423212397157,0.957642744668469];
vbq_def.rfcorr.v2k.RFCorr = true;
% Unknwon protocol
vbq_def.rfcorr.Unknown.tag = 'Unknown protocol. No spoiling correction defined.';
vbq_def.rfcorr.Unknown.RFCorr = false;

%% B1 mapping processing parameters

% For *each* site, the labels corresponding to the available B1 mapping
% protocols must be specified so they are listed as available choices in
% the batch GUI (see tbx_cfg_vbq_crm). Do NOT forget the 'crc'/'fil'/'lren'
% field in the structure. :-)
% The parameters of each of these B1map protocol should be specified in
% their specific substructure, using the protocal name!
% NB: the first label in the list is the default one.
%
% NOTE: 
% all the protocal names MUST 
% - start with a letter, 
% - include only letters, numbers or underscores, i.e. NO space.
% as these names are used to define a structure fieldname with the protocol 
% parameters.

% List B1 protocols available at the CRC
% --------------------------------------
vbq_def.crc.b1_type.labels  = {
    'i3D_EPI_v2b_prisma_crc' % added the 'i' before the '3' to start with a letter...
    'i3D_EPI_v3a_allegra_crc'
    'i3D_EPI_v4a_allegra_crc'
    'i3D_AFI_v4b_n3_allegra_crc'
    'i3D_AFI_v4b_n5_allegra_crc'
    'pre_processed_B1'
    'no_B1_provided'
    }';
vbq_def.crc.b1_type.val  = vbq_def.crc.b1_type.labels(1);
% List B1 protocols available at the FIL
% --------------------------------------
vbq_def.fil.b1_type.labels = {
    'i3D_EPI_v2b_long'
    'i3D_EPI_rect700'
    'pre_processed_B1'
    'no_B1_provided'
    };
vbq_def.fil.b1_type.val = vbq_def.fil.b1_type.labels(1);
% List B1 protocols available at the LREN
% ---------------------------------------
vbq_def.lren.b1_type.labels = {
    'i3D_EPI_v2b_long'
    'i3D_EPI_rect700'
    'pre_processed_B1'
    'no_B1_provided'
    };
vbq_def.lren.b1_type.val = vbq_def.lren.b1_type.labels(1);

% B1 map protocol parameters
% --------------------------
% Each of the B1map protocols, for *all* the sites, are defined in a
% separate substructure.
%
% 1) 'i3D_EPI_v3a_allegra_crc'
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.data    = 'EPI'; 
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.avail   = true; 
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.procreq = true; 
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.T1 = 1192; % ms, strictly valid only at 3T
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.eps = 0.0001;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.beta = 115:-5:65;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.TM = 45;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.Nonominalvalues = 5;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.EchoSpacing = 540e-3;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.nPEacq = 36; % [36] if 75% FoV & PF phase; [48] if only one of these
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.PEDIR = 2; % [2] for R>>L; [1] for A>>P
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b1proc.blipDIR = 1; % +1 for R>>L; -1 for A>>P
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.shorTE = 4.92; % ms
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.longTE = 7.38; % ms
% due to higher level of field inhomogeneity on the allegra,
% the threshold must be set higher otherwise no B1 map for area
% around the brain stem. Trio settings: 110Hz - Allegra
% settings: 300 Hz (ebalteau 20120619)
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.HZTHRESH = 300;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.ERODEB1 = 1;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.PADB1 = 3 ;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.match_vdm = 1;
vbq_def.b1map.i3D_EPI_v3a_allegra_crc.b0proc.b0maskbrain = 0;
% 2) 'i3D_EPI_v4a_allegra_crc'
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.data    = 'EPI'; 
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.avail   = true; 
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.procreq = true; 
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.T1 = 1192; % ms, strictly valid only at 3T
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.eps = 0.0001;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.beta = 140:-7.5:65;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.TM = 40;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.Nonominalvalues = 5;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.EchoSpacing = 330e-3;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.nPEacq = 48; % [36] if 75% FoV & PF phase; [48] if only one of these
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.PEDIR = 2; % [2] for R>>L; [1] for A>>P
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b1proc.blipDIR = 1; % +1 for R>>L; -1 for A>>P
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.shorTE = 4.92; % ms
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.longTE = 7.38; % ms
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.HZTHRESH = 300;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.ERODEB1 = 1;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.PADB1 = 3 ;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.match_vdm = 1;
vbq_def.b1map.i3D_EPI_v4a_allegra_crc.b0proc.b0maskbrain = 0;
% 3) 'i3D_EPI_v2b_prisma_crc' % used on the Prisma (Lausanne & Liege)
% NB: in Liege, maskbrain is set to zero because segmentation
% used in masking procedure diverges for unknown reason...
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.data    = 'EPI'; 
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.avail   = true; 
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.procreq = true; 
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.T1 = 1192; % ms, strictly valid only at 3T
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.eps = 0.0001;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.beta = 115:-5:65;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.TM = 31.2;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.Nonominalvalues = 5;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.EchoSpacing = 540e-3;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.nPEacq = 24; % [36] if 75% FoV & PF phase; [48] if only one of these
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.PEDIR = 2; % [2] for R>>L; [1] for A>>P
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b1proc.blipDIR = 1; % +1 for R>>L; -1 for A>>P
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.shorTE = 10; % ms
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.longTE = 12.46; % ms
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.HZTHRESH = 110;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.ERODEB1 = 1;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.PADB1 = 3 ;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.match_vdm = 1;
vbq_def.b1map.i3D_EPI_v2b_prisma_crc.b0proc.b0maskbrain = 0;
% 4) 'i3D_EPI_v2b_long'
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.data    = 'EPI'; 
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.avail   = true; 
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.procreq = true; 
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.T1 = 1192; % ms, strictly valid only at 3T
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.eps = 0.0001;
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.beta = 135:-5:65;
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.TM = 33.24;
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.Nonominalvalues = 5;
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.EchoSpacing = 540e-3;
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.nPEacq = 24; % [36] if 75% FoV & PF phase; [48] if only one of these
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.PEDIR = 2; % [2] for R>>L; [1] for A>>P
vbq_def.b1map.i3D_EPI_v2b_long.b1proc.blipDIR = 1; % +1 for R>>L; -1 for A>>P
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.shorTE = 10; % ms
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.longTE = 12.46; % ms
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.HZTHRESH = 110;
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.ERODEB1 = 1;
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.PADB1 = 3 ;
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.match_vdm = 1;
vbq_def.b1map.i3D_EPI_v2b_long.b0proc.b0maskbrain = 1;
% 5) 'i3D_EPI_rect700'
vbq_def.b1map.i3D_EPI_rect700.b1proc.data    = 'EPI'; 
vbq_def.b1map.i3D_EPI_rect700.b1proc.procreq = true; 
vbq_def.b1map.i3D_EPI_rect700.b1proc.procreq = true; 
vbq_def.b1map.i3D_EPI_rect700.b1proc.T1 = 1192; % ms, strictly valid only at 3T
vbq_def.b1map.i3D_EPI_rect700.b1proc.eps = 0.0001;
vbq_def.b1map.i3D_EPI_rect700.b1proc.beta = 80:5:100;
vbq_def.b1map.i3D_EPI_rect700.b1proc.TM = 33.53;
vbq_def.b1map.i3D_EPI_rect700.b1proc.Nonominalvalues = 5;
vbq_def.b1map.i3D_EPI_rect700.b1proc.EchoSpacing = 540e-3;
vbq_def.b1map.i3D_EPI_rect700.b1proc.nPEacq = 24; % [36] if 75% FoV & PF phase; [48] if only one of these
vbq_def.b1map.i3D_EPI_rect700.b1proc.PEDIR = 2; % [2] for R>>L; [1] for A>>P
vbq_def.b1map.i3D_EPI_rect700.b1proc.blipDIR = 1; % +1 for R>>L; -1 for A>>P
vbq_def.b1map.i3D_EPI_rect700.b0proc.shorTE = 10; % ms
vbq_def.b1map.i3D_EPI_rect700.b0proc.longTE = 12.46; % ms
vbq_def.b1map.i3D_EPI_rect700.b0proc.HZTHRESH = 110;
vbq_def.b1map.i3D_EPI_rect700.b0proc.ERODEB1 = 1;
vbq_def.b1map.i3D_EPI_rect700.b0proc.PADB1 = 3 ;
vbq_def.b1map.i3D_EPI_rect700.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
vbq_def.b1map.i3D_EPI_rect700.b0proc.match_vdm = 1;
vbq_def.b1map.i3D_EPI_rect700.b0proc.b0maskbrain = 1;
% 6) 'i3D_AFI_v4b_n5_allegra_crc'
vbq_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.data    = 'AFI'; 
vbq_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.avail   = true; 
vbq_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.procreq = true; 
vbq_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.TR2TR1ratio = 5;
vbq_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.alphanom = 60;
% 7) 'i3D_AFI_v4b_n3_allegra_crc'
vbq_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.data    = 'AFI'; 
vbq_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.avail   = true; 
vbq_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.procreq = true; 
vbq_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.TR2TR1ratio = 3;
vbq_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.alphanom = 60;
% 8) 'pre_processed_B1'
vbq_def.b1map.pre_processed_B1.b1proc.avail   = true;
vbq_def.b1map.pre_processed_B1.b1proc.procreq = false;
%9) 'no_B1_provided'
vbq_def.b1map.no_B1_provided.b1proc.procreq = false;
vbq_def.b1map.no_B1_provided.b1proc.avail   = false;

end

%% %%%%%%%%%%%%%%%%% Centre specific parameters %%%%%%%%%%%%%%%%%%%%%%%
%
% Note the centre specific defaults structure
% - they should all have the *same* organization, otherwise crashes could
%    occure when trying to access the default value in the batch.
% - they should NOT have anything similar with the 'global' defaults.
%   For example do NOT define "vbq_def.TE = 20" and "vbq_def.fil.TE = 30"
%   fields as the latter would NEVER be used

%% Specific parameters, CRC
% Examples:
% vbq_def.crc.TR  = 3;  % in sec
% vbq_def.crc.TE1 = 50; % in ms
% vbq_def.crc.TE2 = 80; % in ms
% vbq_def.crc.cset1.val1 = 12; % in ms
% vbq_def.crc.cset1.val2 = 34; % in ms
% vbq_def.crc.cset2.val1 = 56; % in ms
% vbq_def.crc.cset2.val2 = 78; % in ms


%% Specific parameters, FIL
% % Examples:
% vbq_def.fil.TR  = 2;  % in sec
% vbq_def.fil.TE1 = 30; % in ms
% vbq_def.fil.TE2 = 60; % in ms
% vbq_def.fil.cset1.val1 = 21; % in ms
% vbq_def.fil.cset1.val2 = 43; % in ms
% vbq_def.fil.cset2.val1 = 65; % in ms
% vbq_def.fil.cset2.val2 = 87; % in ms


%% Specific parameters, LReN
% % Examples:
% vbq_def.lren.TR  = 2.5;  % in sec
% vbq_def.lren.TE1 = 40; % in ms
% vbq_def.lren.TE2 = 70; % in ms
% vbq_def.lren.cset1.val1 = 13; % in ms
% vbq_def.lren.cset1.val2 = 24; % in ms
% vbq_def.lren.cset2.val1 = 47; % in ms
% vbq_def.lren.cset2.val2 = 68; % in ms


