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

% Note: 
% it's advised to organize the defaults by "families", i.e. substructures.

%% 
global vbq_def

%% %%%%%%%%%%%%%%%%%%%%% Global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specifying the lab
vbq_def.centre = 'crc' ; % can be 'fil' or 'lren' or 'crc'

% %% Some other global parameters
% vbq_def.param1 = '123' ;
% vbq_def.param2 = '456' ;
% 
% %% Specific set of parameters, SET1
% vbq_def.set1.prefix = 'a'; % when a char is specified
% vbq_def.set1.val1   = 1; % setting a single value
% 
% %% Specific set of parameters, SET2
% vbq_def.set2.prefix = 'b'; % when a char is specified
% vbq_def.set2.val1   = 2; % setting a single value


%% %%%%%%%%%%%%%%%%% Common processing parameters %%%%%%%%%%%%%%%%%%%%%
% These parameters are either parameters that are fixed for all sites or
% recommended values. They can also be changed in a site-specific way at
% run-time. 

%% Processing of PD maps
vbq_def.PDproc.PDmap = 1;       %Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
vbq_def.PDproc.WBMaskTh = 0.1;  %Threshold for calculation of whole-brain mask from TPMs
vbq_def.PDproc.WMMaskTh = 0.95; %Threshold for calculation of white-matter mask from TPMs
vbq_def.PDproc.biasreg = 10^(-5);
vbq_def.PDproc.biasfwhm = 50;

%% Threshold values for saving of the qMRI maps
vbq_def.qMRI_maps_thresh.R1 = 2000;
vbq_def.qMRI_maps_thresh.A = 10^5;
vbq_def.qMRI_maps_thresh.R2s = 10;
vbq_def.qMRI_maps_thresh.MTR = 50;
vbq_def.qMRI_maps_thresh.MTR_synt = 50;
vbq_def.qMRI_maps_thresh.MT = 5;

%% MPM acquisition parameters and RF spoiling correction parameters
% these value are initialised with defaults (v2k protocol - Prisma) for the
% first pass through this script only. They're updated at run-time with
% actual acquisition values (see vbq_mpr_b0_b1.m). The coefficients for
% R.Deichmann steady state correction are also determined and stored as
% part of the defaults parameters for the current processing so they can be
% stored for reccord. 
if ~isfield(vbq_def,'MPMacq')
    vbq_def.MPMacq.TE_mtw = 2.34;
    vbq_def.MPMacq.TE_t1w = 2.34;
    vbq_def.MPMacq.TE_pdw = 2.34;
    vbq_def.MPMacq.TR_mtw = 24.5;
    vbq_def.MPMacq.TR_t1w = 24.5;
    vbq_def.MPMacq.TR_pdw = 24.5;
    vbq_def.MPMacq.fa_mtw = 6;
    vbq_def.MPMacq.fa_t1w = 21;
    vbq_def.MPMacq.fa_pdw = 6;
end

% Antoine Lutti 15/01/09
% Correction parameters used in vbq_MTProt to correct for imperfect RF 
% spoiling when a B1 map is loaded. Correction based on Preibisch and 
% Deichmann's paper MRM 61:125-135 (2009). The values for P2_a and P2_b 
% below were obtained using the code supplied by R. Deichmann with the 
% experimentalparameters used to get our PDw and T1w images. Correction 
% parameters were calculated for the following parameter sets using 
% T2 = 64 ms at 3T.

if ((vbq_def.MPMacq.TR_pdw == 23.7) && ...
        (vbq_def.MPMacq.TR_t1w == 18.7) && ...
        (vbq_def.MPMacq.fa_pdw == 6) && ...
        (vbq_def.MPMacq.fa_t1w == 20))
    % 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
    vbq_def.rfcorr.tag = 'Classic FIL protocol';
    vbq_def.rfcorr.P2_a = [78.9228195006542,-101.113338489192,47.8783287525126];
    vbq_def.rfcorr.P2_b = [-0.147476233142129,0.126487385091045,0.956824374979504];
    vbq_def.rfcorr.RFCorr = true;
elseif ((vbq_def.MPMacq.TR_pdw == 24.5) && ...
        (vbq_def.MPMacq.TR_t1w == 24.5) && ...
        (vbq_def.MPMacq.fa_pdw == 5) && ...
        (vbq_def.MPMacq.fa_t1w == 29))
    % 2) new FIL/Helms protocol
    % PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
    vbq_def.rfcorr.tag = 'New FIL/Helms protocol';
    vbq_def.rfcorr.P2_a = [93.455034845930480,-120.5752858196904,55.911077913369060];
    vbq_def.rfcorr.P2_b = [-0.167301931434861,0.113507432776106,0.961765216743606];
    vbq_def.rfcorr.RFCorr = true;
elseif ((vbq_def.MPMacq.TR_pdw == 24.0) && ...
        (vbq_def.MPMacq.TR_t1w == 19.0) && ...
        (vbq_def.MPMacq.fa_pdw == 6) && ...
        (vbq_def.MPMacq.fa_t1w == 20))
    % 3) Siemens product sequence protocol used in Lausanne (G Krueger)
    %PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
    vbq_def.rfcorr.tag = 'Siemens product Lausanne (GK) protocol';
    vbq_def.rfcorr.P2_a = [67.023102027100880,-86.834117103841540,43.815818592349870];
    vbq_def.rfcorr.P2_b = [-0.130876849571103,0.117721807209409,0.959180058389875];
    vbq_def.rfcorr.RFCorr = true;
elseif ((vbq_def.MPMacq.TR_pdw ==  23.7) && ...
        (vbq_def.MPMacq.TR_t1w == 23.7) && ...
        (vbq_def.MPMacq.fa_pdw == 6) && ...
        (vbq_def.MPMacq.fa_t1w == 28))
    % 4) High-res (0.8mm) FIL protocol:
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
    vbq_def.rfcorr.tag = 'High-res FIL protocol';
    vbq_def.rfcorr.P2_a = [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
    vbq_def.rfcorr.P2_b = [-0.218804328507184,0.178745853134922,0.939514554747592];
    vbq_def.rfcorr.RFCorr = true;
elseif ((vbq_def.MPMacq.TR_pdw == 25.25) && ...
        (vbq_def.MPMacq.TR_t1w == 25.25) && ...
        (vbq_def.MPMacq.fa_pdw == 5) && ...
        (vbq_def.MPMacq.fa_t1w == 29))
    % 4)NEW  High-res (0.8mm) FIL protocol:
    % PD-weighted: TR=25.25ms; a=5deg; T1-weighted: TR=TR=25.25ms; a=29deg
    vbq_def.rfcorr.tag = 'High-res FIL protocol';
    vbq_def.rfcorr.P2_a = [88.8623036106612,-114.526218941363,53.8168602253166];
    vbq_def.rfcorr.P2_b = [-0.132904017579521,0.113959390779008,0.960799295622202];
    vbq_def.rfcorr.RFCorr = true;
elseif ((vbq_def.MPMacq.TR_pdw == 24.5) && ...
        (vbq_def.MPMacq.TR_t1w == 24.5) && ...
        (vbq_def.MPMacq.fa_pdw == 6) && ...
        (vbq_def.MPMacq.fa_t1w == 21))
    % 5)NEW  1mm protocol - seq version v2k:
    % PD-weighted: TR=24.5ms; a=6deg; T1-weighted: TR=24.5ms; a=21deg
    vbq_def.rfcorr.tag = 'v2k protocol';
    vbq_def.rfcorr.P2_a = [71.2817617982844,-92.2992876164017,45.8278193851731];
    vbq_def.rfcorr.P2_b = [-0.137859046784839,0.122423212397157,0.957642744668469];
    vbq_def.rfcorr.RFCorr = true;
else
    warning('Warning!!! Spoiling correction not defined for this protocol. No correction being applied.');
    vbq_def.rfcorr.tag = 'Unknown protocol. No spoiling correction defined.';
    vbq_def.rfcorr.RFCorr = false;
end

%% B1 mapping processing parameters
% Although part of the parameters are site-specific and others are common
% to all sites, a single structure is defined. Site-specific values must
% overwrite default, common values.

% +++++++ COMMON AND DEFAULT PARAMETERS +++++++
vbq_def.b1map.T1 = 1192; %ms, strictly valid only at 3T
vbq_def.b1map.eps = 0.0001;
vbq_def.b1map.b1_type.labels = {'3D_EPI_v2b'};
vbq_def.b1map.avail = 1; % by default assumes b1 mapping data are available
vbq_def.b1map.procreq = 1; % ... and processing required to create b1 map.

% +++++++ SITE SPECIFIC PARAMETERS +++++++
% The following parameters are site-specific but also depending on the
% b1_type. Therefore they'll be assigned at runtime according to the exact
% b1_type.val (initially undefined) by calling:
%       vbq_get_defaults('b1map.b1_type.val',job.b1_type); 
% (to specify the b1_type)
% and:
%       vbq_defaults; 
% (to force running through the code below with new b1_type)
% from run_vbq_b1map.
%
% For each site, the labels corresponding to the available B1 mapping
% protocols must be specified so they are listed as available choices in
% the batch GUI (see tbx_cfg_vbq). NB: the first label in the list is the
% default type.
switch vbq_def.centre
    case 'crc'
        vbq_def.b1map.b1_type.labels  = {
            '3D_EPI_v2b_prisma_crc'
            '3D_EPI_v3a_allegra_crc'
            '3D_EPI_v4a_allegra_crc'
            '3D_AFI_v4b_n3_allegra_crc'
            '3D_AFI_v4b_n5_allegra_crc'
            'pre_processed_B1'
            'no_B1_provided'
            }';
    case 'fil'
        
    case 'lren'
        
end

% The specific B1 processing values are assigned as soon as the b1_type is
% known (which might not be the case at the first call of the vbq_defaults
% function): 
if isfield(vbq_def.b1map.b1_type,'val')
    switch vbq_def.b1map.b1_type.val
        case '3D_EPI_v3a_allegra_crc'
            vbq_def.b1map.beta = 115:-5:65;
            vbq_def.b1map.TM = 45;
            vbq_def.b1map.Nonominalvalues = 5;
            vbq_def.b1map.EchoSpacing = 540e-3;
            vbq_def.b1map.nPEacq = 36; % [36] if 75% FoV & PF phase; [48] if only one of these
            vbq_def.b1map.PEDIR = 2; % [2] for R>>L; [1] for A>>P
            vbq_def.b1map.blipDIR = 1; % +1 for R>>L; -1 for A>>P
            vbq_def.b0proc.shorTE = 4.92; % ms
            vbq_def.b0proc.longTE = 7.38; % ms
            % due to higher level of field inhomogeneity on the allegra,
            % the threshold must be set higher otherwise no B1 map for area
            % around the brain stem. Trio settings: 110Hz - Allegra
            % settings: 300 Hz (ebalteau 20120619)
            vbq_def.b0proc.HZTHRESH = 300;
            vbq_def.b0proc.ERODEB1 = 1;
            vbq_def.b0proc.PADB1 = 3 ;
            vbq_def.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
            vbq_def.b0proc.match_vdm = 1;
            vbq_def.b0proc.b0maskbrain = 0;
        case '3D_EPI_v4a_allegra_crc'
            vbq_def.b1map.beta = 140:-7.5:65;
            vbq_def.b1map.TM = 40;
            vbq_def.b1map.Nonominalvalues = 5;
            vbq_def.b1map.EchoSpacing = 330e-3;
            vbq_def.b1map.nPEacq = 48; % [36] if 75% FoV & PF phase; [48] if only one of these
            vbq_def.b1map.PEDIR = 2; % [2] for R>>L; [1] for A>>P
            vbq_def.b1map.blipDIR = 1; % +1 for R>>L; -1 for A>>P
            vbq_def.b0proc.shorTE = 4.92; % ms
            vbq_def.b0proc.longTE = 7.38; % ms
            vbq_def.b0proc.HZTHRESH = 300;
            vbq_def.b0proc.ERODEB1 = 1;
            vbq_def.b0proc.PADB1 = 3 ;
            vbq_def.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
            vbq_def.b0proc.match_vdm = 1;
            vbq_def.b0proc.b0maskbrain = 0;
        case '3D_EPI_v2b_prisma_crc' % used on the Prisma (Lausanne & Liege)
            % NB: in Liege, maskbrain is set to zero because segmentation
            % used in masking procedure diverges for unknown reason...
            vbq_def.b1map.beta = 115:-5:65;
            vbq_def.b1map.TM = 31.2;
            vbq_def.b1map.Nonominalvalues = 5;
            vbq_def.b1map.EchoSpacing = 540e-3;
            vbq_def.b1map.nPEacq = 24; % [36] if 75% FoV & PF phase; [48] if only one of these
            vbq_def.b1map.PEDIR = 2; % [2] for R>>L; [1] for A>>P
            vbq_def.b1map.blipDIR = 1; % +1 for R>>L; -1 for A>>P
            vbq_def.b0proc.shorTE = 10; % ms
            vbq_def.b0proc.longTE = 12.46; % ms        
            vbq_def.b0proc.HZTHRESH = 110;
            vbq_def.b0proc.ERODEB1 = 1;
            vbq_def.b0proc.PADB1 = 3 ;
            vbq_def.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
            vbq_def.b0proc.match_vdm = 1;
            vbq_def.b0proc.b0maskbrain = 0;
        case '3D_EPI_v2b_long'
            vbq_def.b1map.beta = 135:-5:65;
            vbq_def.b1map.TM = 33.24;
            vbq_def.b1map.Nonominalvalues = 5;
            vbq_def.b1map.EchoSpacing = 540e-3;
            vbq_def.b1map.nPEacq = 24; % [36] if 75% FoV & PF phase; [48] if only one of these
            vbq_def.b1map.PEDIR = 2; % [2] for R>>L; [1] for A>>P
            vbq_def.b1map.blipDIR = 1; % +1 for R>>L; -1 for A>>P
            vbq_def.b0proc.shorTE = 10; % ms
            vbq_def.b0proc.longTE = 12.46; % ms        
            vbq_def.b0proc.HZTHRESH = 110;
            vbq_def.b0proc.ERODEB1 = 1;
            vbq_def.b0proc.PADB1 = 3 ;
            vbq_def.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
            vbq_def.b0proc.match_vdm = 1;
            vbq_def.b0proc.b0maskbrain = 1;
        case '3D_EPI_rect700'
            vbq_def.b1map.beta = 80:5:100;
            vbq_def.b1map.TM = 33.53;
            vbq_def.b1map.Nonominalvalues = 5;
            vbq_def.b1map.EchoSpacing = 540e-3;
            vbq_def.b1map.nPEacq = 24; % [36] if 75% FoV & PF phase; [48] if only one of these
            vbq_def.b1map.PEDIR = 2; % [2] for R>>L; [1] for A>>P
            vbq_def.b1map.blipDIR = 1; % +1 for R>>L; -1 for A>>P
            vbq_def.b0proc.shorTE = 10; % ms
            vbq_def.b0proc.longTE = 12.46; % ms        
            vbq_def.b0proc.HZTHRESH = 110;
            vbq_def.b0proc.ERODEB1 = 1;
            vbq_def.b0proc.PADB1 = 3 ;
            vbq_def.b0proc.B1FWHM = 8; %For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
            vbq_def.b0proc.match_vdm = 1;
            vbq_def.b0proc.b0maskbrain = 1;
        case '3D_AFI_v4b_n5_allegra_crc'
            vbq_def.b1map.TR2TR1ratio = 5;
            vbq_def.b1map.alphanom = 60;
        case '3D_AFI_v4b_n3_allegra_crc'
            vbq_def.b1map.TR2TR1ratio = 3;
            vbq_def.b1map.alphanom = 60;
        case 'pre_processed_B1'
            vbq_def.b1map.procreq = 0;
        case 'no_B1_provided'
            vbq_def.b1map.procreq = 0;
            vbq_def.b1map.avail = 0;
        otherwise % 'no_B1_provided'
            vbq_def.b1map.procreq = 0;
            vbq_def.b1map.avail = 0;
            warning('Warning!!! Unknown type of B1 map data. No B1 correction will be performed');
    end
end

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


