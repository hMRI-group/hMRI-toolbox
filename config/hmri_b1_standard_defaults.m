function hmri_b1_standard_defaults
% Sets the defaults for B1 bias correction, part of the hMRI toolbox.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% DON'T MODIFY THIS FILE, IT CONTAINS B1 PROCESSING DEFAULTS PARAMETERS.
% Please refer to hMRI-Toolbox\config\local\hmri_b1_local_defaults.m to
% customise B1 processing parameters.  
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% FORMAT hmri_b1_standard_defaults
%__________________________________________________________________________
%
% PURPOSE
% To set the reference defaults parameters to process B1 data for B1 bias
% correction (multiparameter mapping protocol of the hMRI-Toolbox).
%
% For B1 map calculation parameters customization (available for 3DEPI,
% 3DAFI and UNICORT methods), use "hmri_b1_local_defaults.m". Make a copy
% of the latter file with meaningful name, modify as desired and select as
% defaults B1 file in the MapCreation Batch GUI.
%
% WARNING
% Modification of the defaults parameters may impair the the integrity of
% the toolbox, leading to unexpected behaviour. Only recommended for expert
% users. 
%
%__________________________________________________________________________
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium
%__________________________________________________________________________

% Global hmri_def variable used across the whole toolbox
global hmri_def

%--------------------------------------------------------------------------
% B1 mapping processing parameters 
%--------------------------------------------------------------------------
% Default parameters are set below for each type of B1 processing.
% For acquisition parameters, default values are a fallback solution for B1
% data processing when no metadata are available. Use of metadata is
% recommended to retrieve site- & protocol-specific parameters and ensure
% appropriate data handling and processing.

% 'i3D_AFI'
hmri_def.b1map.i3D_AFI.b1type = 'i3D_AFI'; 
hmri_def.b1map.i3D_AFI.b1avail = true; 
hmri_def.b1map.i3D_AFI.procreq = true; 
hmri_def.b1map.i3D_AFI.b1acq.TR2TR1ratio = 5;
hmri_def.b1map.i3D_AFI.b1acq.alphanom = 60;

% 'pre_processed_B1'
hmri_def.b1map.pre_processed_B1.b1type = 'pre_processed_B1'; 
hmri_def.b1map.pre_processed_B1.b1avail   = true;
hmri_def.b1map.pre_processed_B1.procreq = false;

% 'no_B1_correction'
hmri_def.b1map.no_B1_correction.b1type = 'no_B1_correction'; 
hmri_def.b1map.no_B1_correction.b1avail   = false;
hmri_def.b1map.no_B1_correction.procreq = false;

% UNICORT
hmri_def.b1map.UNICORT.b1type = 'UNICORT'; 
hmri_def.b1map.UNICORT.procreq = true;
hmri_def.b1map.UNICORT.b1avail   = false;
hmri_def.b1map.UNICORT.procpar.reg = 10^-3;
hmri_def.b1map.UNICORT.procpar.FWHM = 60;
hmri_def.b1map.UNICORT.procpar.thr = 5;

% 'i3D_EPI'
hmri_def.b1map.i3D_EPI.b1type = 'i3D_EPI'; 
hmri_def.b1map.i3D_EPI.b1avail   = true; 
hmri_def.b1map.i3D_EPI.procreq = true; 
% b0&b1-processing
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
% b1-acquisition
% hmri_def.b1map.i3D_EPI.b1acq.beta = 115:-5:65;
% hmri_def.b1map.i3D_EPI.b1acq.TM = 31.2;
% for long acqu
hmri_def.b1map.i3D_EPI.b1acq.beta = 135:-5:65;
hmri_def.b1map.i3D_EPI.b1acq.TM = 33.24;
hmri_def.b1map.i3D_EPI.b1acq.tert = 540e-3*24; % EchoSpacing * numberPElines
hmri_def.b1map.i3D_EPI.b1acq.blipDIR = 1;
% b0-acquisition
hmri_def.b1map.i3D_EPI.b0acq.shortTE = 10; % ms
hmri_def.b1map.i3D_EPI.b0acq.longTE = 12.46; % ms
hmri_def.b1map.i3D_EPI.b0acq.iformat = 'PM'; % ms

% 'tfl_b1_map'
hmri_def.b1map.tfl_b1_map.b1type = 'tfl_b1_map'; 
hmri_def.b1map.tfl_b1_map.b1avail   = true; 
hmri_def.b1map.tfl_b1_map.procreq = true; 

% 'rf_map'
hmri_def.b1map.rf_map.b1type = 'rf_map'; 
hmri_def.b1map.rf_map.b1avail   = true; 
hmri_def.b1map.rf_map.procreq = true; 

end