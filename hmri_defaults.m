function hmri_defaults
% Sets the defaults which are used by the hMRI toolbox.
%
% FORMAT hmri_defaults
%_______________________________________________________________________
%
% This file can be customised to any the site/person own setup.
% Individual users can make copies which can be stored on their own
% matlab path. Make sure then that your 'hmri_defaults' is the first one
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

%%
global hmri_def

%% %%%%%%%%%%%%%%%%%%%%% Global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specifying the lab
hmri_def.centre = 'cbs' ; % can be 'fil', 'lren', 'crc', 'sciz' or 'cbs'

%% %%%%%%%%%%%%%%%%% Common processing parameters %%%%%%%%%%%%%%%%%%%%%
% These parameters are either parameters that are fixed for all sites or
% recommended values. They can also be changed in a site-specific way at
% run-time.

hmri_def.R2sOLS = 1; % Create an Ordinary Least Squares R2* map?
hmri_def.json = struct('extended',false,'separate',true,'anonym','none',...
    'overwrite',true); % settings for JSON metadata

%% Processing of PD maps
hmri_def.PDproc.PDmap    = 1;    % Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
hmri_def.PDproc.WBMaskTh = 0.1;  % Threshold for calculation of whole-brain mask from TPMs
hmri_def.PDproc.WMMaskTh = 0.95; % Threshold for calculation of white-matter mask from TPMs
hmri_def.PDproc.biasreg  = 10^(-5);
hmri_def.PDproc.biasfwhm = 50;
hmri_def.PDproc.nr_echoes_forA =1;

%% UNICORT processing
hmri_def.unicort.reg = 10^-3;
hmri_def.unicort.FWHM = 60;
hmri_def.unicort.thr = 2; % TL: 2 for sciz & cbs with SIEMENS 3T Skyra fit
                          % otherwise: 5

hmri_def.qMRI_maps.QA          = 1;%creates a matlab structure containing markers of data quality
hmri_def.qMRI_maps.ACPCrealign = 1;%realigns qMRI maps to MNI

%% Threshold values for saving of the qMRI maps
hmri_def.qMRI_maps_thresh.R1       = 2000;
hmri_def.qMRI_maps_thresh.A        = 10^5;
hmri_def.qMRI_maps_thresh.R2s      = 10;
hmri_def.qMRI_maps_thresh.MTR      = 50;
hmri_def.qMRI_maps_thresh.MTR_synt = 50;
hmri_def.qMRI_maps_thresh.MT       = 15; % TL: 15 for cbs & sciz with SIEMENS 3T Skyra
                                         % original: 5 

%% MPM acquisition parameters and RF spoiling correction parameters
% these value are initialised with defaults (v2k protocol - Prisma) for the
% first pass through this script only. They're updated at run-time with
% actual acquisition values (see hmri_MTProt.m).
% The coefficients for R.Deichmann steady state correction are also
% determined and stored as part of the defaults parameters for the current
% processing so they can be stored for reccord.

% Default MPMacq values to begin with
%---------------------------------------------------------------------
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
%---------------------------------------------------------------------
% Antoine Lutti 15/01/09
% Correction parameters used in hmri_MTProt to correct for imperfect RF
% spoiling when a B1 map is loaded. Correction based on Preibisch and
% Deichmann's paper MRM 61:125-135 (2009). The values for P2_a and P2_b
% below were obtained using the code supplied by R. Deichmann with the
% experimentalparameters used to get our PDw and T1w images. Correction
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

%% B1 mapping processing parameters

% For *each* site, the labels corresponding to the available B1 mapping
% protocols must be specified so they are listed as available choices in
% the batch GUI (see tbx_cfg_hmri_crm). Do NOT forget the 'crc'/'fil'/'lren'
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

% List B1 protocols available at the CBS
% --------------------------------------
hmri_def.cbs.b1_type.labels  = {
    'i3D_EPI'
    'tfl_b1_map'
    'rf_map'
    'no_B1_correction'
    'pre_processed_B1'
    'UNICORT'
    }';
hmri_def.cbs.b1_type.val  = hmri_def.cbs.b1_type.labels(1);

% List B1 protocols available at the CRC
% --------------------------------------
hmri_def.crc.b1_type.labels  = {
    'i3D_EPI' % added the 'i' before the '3' to start with a letter...
    'i3D_AFI_v4b_n3_allegra_crc'
    'i3D_AFI_v4b_n5_allegra_crc'
    'no_B1_correction'
    'pre_processed_B1'
    'UNICORT'
    }';
hmri_def.crc.b1_type.val  = hmri_def.crc.b1_type.labels(1);

% List B1 protocols available at the FIL
% --------------------------------------
hmri_def.fil.b1_type.labels = {
    'i3D_EPI'
    'no_B1_correction'
    'pre_processed_B1'
    'UNICORT'
    };
hmri_def.fil.b1_type.val = hmri_def.fil.b1_type.labels(1);

% List B1 protocols available at the LREN
% ---------------------------------------
hmri_def.lren.b1_type.labels = {
    'i3D_EPI'
    'no_B1_correction'
    'pre_processed_B1'
    'UNICORT'
    };
hmri_def.lren.b1_type.val = hmri_def.lren.b1_type.labels(1);

% List B1 protocols available at the Spinal Cord Injury Lab, Zurich
% -----------------------------------------------------------------
hmri_def.sciz.b1_type.labels  = {
    'tfl_b1map'
    'rf_map'
    'no_B1_correction'
    'pre_processed_B1'
    'UNICORT'
    }';
hmri_def.sciz.b1_type.val  = hmri_def.sciz.b1_type.labels(1);

% B1 map protocol parameters
% --------------------------
% Each of the B1map protocols, for *all* the sites, are defined in a
% separate substructure.
%
% 'i3D_AFI_v4b_n5_allegra_crc'
hmri_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.data    = 'AFI'; 
hmri_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.avail   = true; 
hmri_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.procreq = true; 
hmri_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.TR2TR1ratio = 5;
hmri_def.b1map.i3D_AFI_v4b_n5_allegra_crc.b1proc.alphanom = 60;

% 'i3D_AFI_v4b_n3_allegra_crc'
hmri_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.data    = 'AFI'; 
hmri_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.avail   = true; 
hmri_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.procreq = true; 
hmri_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.TR2TR1ratio = 3;
hmri_def.b1map.i3D_AFI_v4b_n3_allegra_crc.b1proc.alphanom = 60;

% 'pre_processed_B1'
hmri_def.b1map.pre_processed_B1.b1proc.avail   = true;
hmri_def.b1map.pre_processed_B1.b1proc.procreq = false;

% 'no_B1_correction'
hmri_def.b1map.no_B1_correction.b1proc.procreq = false;
hmri_def.b1map.no_B1_correction.b1proc.avail   = false;

% UNICORT
hmri_def.b1map.UNICORT.b1proc.procreq = true;
hmri_def.b1map.UNICORT.b1proc.avail   = false;

% 'i3D_EPI'
hmri_def.b1map.i3D_EPI.data    = 'EPI'; 
hmri_def.b1map.i3D_EPI.avail   = true; 
hmri_def.b1map.i3D_EPI.procreq = true; 

hmri_def.b1map.i3D_EPI.b1proc.T1 = 1192; % ms, strictly valid only at 3T
hmri_def.b1map.i3D_EPI.b1proc.eps = 0.0001;
hmri_def.b1map.i3D_EPI.b1proc.Nonominalvalues = 5;
hmri_def.b1map.i3D_EPI.b1proc.HZTHRESH = 110;
hmri_def.b1map.i3D_EPI.b1proc.SDTHRESH = 5;
hmri_def.b1map.i3D_EPI.b1proc.ERODEB1 = 1;
hmri_def.b1map.i3D_EPI.b1proc.PADB1 = 3 ;
hmri_def.b1map.i3D_EPI.b1proc.B1FWHM = 8; % For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
hmri_def.b1map.i3D_EPI.b1proc.match_vdm = 1;
hmri_def.b1map.i3D_EPI.b1proc.b0maskbrain = 1;

hmri_def.b1map.i3D_EPI.b1acq.beta = 115:-5:65;
hmri_def.b1map.i3D_EPI.b1acq.TM = 31.2;
hmri_def.b1map.i3D_EPI.b1acq.EchoSpacing = 540e-3;
hmri_def.b1map.i3D_EPI.b1acq.nPEacq = 24;
hmri_def.b1map.i3D_EPI.b1acq.blipDIR = 1;

hmri_def.b1map.i3D_EPI.b0acq.shortTE = 10; % ms
hmri_def.b1map.i3D_EPI.b0acq.longTE = 12.46; % ms

% 13) 'tfl_b1_map'
hmri_def.b1map.tfl_b1_map.data    = 'TFL'; 
hmri_def.b1map.tfl_b1_map.avail   = true; 
hmri_def.b1map.tfl_b1_map.procreq = true; 
% 14) 'rf_map'
hmri_def.b1map.rf_map.data    = 'RFmap'; 
hmri_def.b1map.rf_map.avail   = true; 
hmri_def.b1map.rf_map.procreq = true; 

end

%% %%%%%%%%%%%%%%%%% Centre specific parameters %%%%%%%%%%%%%%%%%%%%%%%
%
% Note the centre specific defaults structure
% - they should all have the *same* organization, otherwise crashes could
%    occure when trying to access the default value in the batch.
% - they should NOT have anything similar with the 'global' defaults.
%   For example do NOT define "hmri_def.TE = 20" and "hmri_def.fil.TE = 30"
%   fields as the latter would NEVER be used

%% Specific parameters, CRC
% Examples:
% hmri_def.crc.TR  = 3;  % in sec
% hmri_def.crc.TE1 = 50; % in ms
% hmri_def.crc.TE2 = 80; % in ms
% hmri_def.crc.cset1.val1 = 12; % in ms
% hmri_def.crc.cset1.val2 = 34; % in ms
% hmri_def.crc.cset2.val1 = 56; % in ms
% hmri_def.crc.cset2.val2 = 78; % in ms


%% Specific parameters, FIL
% % Examples:
% hmri_def.fil.TR  = 2;  % in sec
% hmri_def.fil.TE1 = 30; % in ms
% hmri_def.fil.TE2 = 60; % in ms
% hmri_def.fil.cset1.val1 = 21; % in ms
% hmri_def.fil.cset1.val2 = 43; % in ms
% hmri_def.fil.cset2.val1 = 65; % in ms
% hmri_def.fil.cset2.val2 = 87; % in ms


%% Specific parameters, LReN
% % Examples:
% hmri_def.lren.TR  = 2.5;  % in sec
% hmri_def.lren.TE1 = 40; % in ms
% hmri_def.lren.TE2 = 70; % in ms
% hmri_def.lren.cset1.val1 = 13; % in ms
% hmri_def.lren.cset1.val2 = 24; % in ms
% hmri_def.lren.cset2.val1 = 47; % in ms
% hmri_def.lren.cset2.val2 = 68; % in ms
