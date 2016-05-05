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
vbq_def.centre = 'fil' ; % but could be 'fil' or 'lren'

%% Some other global parameters
vbq_def.param1 = '123' ;
vbq_def.param2 = '456' ;

%% Specific set of parameters, SET1
vbq_def.set1.prefix = 'a'; % when a char is specified
vbq_def.set1.val1   = 1; % setting a single value

%% Specific set of parameters, SET2
vbq_def.set2.prefix = 'b'; % when a char is specified
vbq_def.set2.val1   = 2; % setting a single value


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
% Proton density mapping
vbq_def.crc.PDmap=1;%Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
vbq_def.crc.WBMaskTh=0.1;%Threshold for calculation of whole-brain mask from TPMs
vbq_def.crc.WMMaskTh=0.95;%Threshold for calculation of white-matter mask from TPMs
vbq_def.crc.biasreg=10^(-5);
vbq_def.crc.biasfwhm=50;
% Threshold values for saving of the qMRI maps
vbq_def.crc.R1thresh=2000;
vbq_def.crc.Athresh=10^5;
vbq_def.crc.R2sthresh=10;
vbq_def.crc.MTRthresh=50;
vbq_def.crc.MTR_syntthresh=50;
vbq_def.crc.MTthresh=5;
vbq_def.crc.R2sOLS = 0; % Create an Ordinary Least Squares R2* map?

%% Specific parameters, FIL
% % Examples:
% vbq_def.fil.TR  = 2;  % in sec
% vbq_def.fil.TE1 = 30; % in ms
% vbq_def.fil.TE2 = 60; % in ms
% vbq_def.fil.cset1.val1 = 21; % in ms
% vbq_def.fil.cset1.val2 = 43; % in ms
% vbq_def.fil.cset2.val1 = 65; % in ms
% vbq_def.fil.cset2.val2 = 87; % in ms
% Proton density mapping
vbq_def.fil.PDmap=1;%Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
vbq_def.fil.WBMaskTh=0.1;%Threshold for calculation of whole-brain mask from TPMs
vbq_def.fil.WMMaskTh=0.95;%Threshold for calculation of white-matter mask from TPMs
vbq_def.fil.biasreg=10^(-5);
vbq_def.fil.biasfwhm=50;
% Threshold values for saving of the qMRI maps
vbq_def.fil.R1thresh=2000;
vbq_def.fil.Athresh=10^5;
vbq_def.fil.R2sthresh=10;
vbq_def.fil.MTRthresh=50;
vbq_def.fil.MTR_syntthresh=50;
vbq_def.fil.MTthresh=5;
vbq_def.fil.R2sOLS = 0; % Create an Ordinary Least Squares R2* map?

%% Specific parameters, LReN
% % Examples:
% vbq_def.lren.TR  = 2.5;  % in sec
% vbq_def.lren.TE1 = 40; % in ms
% vbq_def.lren.TE2 = 70; % in ms
% vbq_def.lren.cset1.val1 = 13; % in ms
% vbq_def.lren.cset1.val2 = 24; % in ms
% vbq_def.lren.cset2.val1 = 47; % in ms
% vbq_def.lren.cset2.val2 = 68; % in ms
% Proton density mapping
vbq_def.lren.PDmap=1;%Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
vbq_def.lren.WBMaskTh=0.1;%Threshold for calculation of whole-brain mask from TPMs
vbq_def.lren.WMMaskTh=0.95;%Threshold for calculation of white-matter mask from TPMs
vbq_def.lren.biasreg=10^(-5);
vbq_def.lren.biasfwhm=50;
% Threshold values for saving of the qMRI maps
vbq_def.lren.R1thresh=2000;
vbq_def.lren.Athresh=10^5;
vbq_def.lren.R2sthresh=10;
vbq_def.lren.MTRthresh=50;
vbq_def.lren.MTR_syntthresh=50;
vbq_def.lren.MTthresh=5;
vbq_def.lren.R2sOLS = 0; % Create an Ordinary Least Squares R2* map?

%% Specific parameters, Spinal Cord Injury Lab, Zeurich: "sciz"
% Proton density mapping
vbq_def.sciz.PDmap=1;%Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
vbq_def.sciz.WBMaskTh=0.1;%Threshold for calculation of whole-brain mask from TPMs
vbq_def.sciz.WMMaskTh=0.95;%Threshold for calculation of white-matter mask from TPMs
vbq_def.sciz.biasreg=10^(-5);
vbq_def.sciz.biasfwhm=50;
% Threshold values for saving of the qMRI maps
vbq_def.sciz.R1thresh=2000;
vbq_def.sciz.Athresh=10^5;
vbq_def.sciz.R2sthresh=10;
vbq_def.sciz.MTRthresh=50;
vbq_def.sciz.MTR_syntthresh=50;
vbq_def.sciz.MTthresh=15;
vbq_def.sciz.R2sOLS = 0; % Create an Ordinary Least Squares R2* map?

return

