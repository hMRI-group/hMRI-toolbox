function hmri_b1_defaults
% Sets the defaults for B1 bias correction, part of the hMRI toolbox.
%
% FORMAT hmri_b1_defaults
%__________________________________________________________________________
%
% This file can be customised to match any site/person own B1 mapping
% protocols.
% Individual users can make copies of this file with appropriate 
% parameters and meaningful name. The customised file can be selected in
% the batch GUI to process B1 bias correction with appropriate parameters.
% Note that this is a fallback solution for situations where acquisition
% and processing parameters cannot be derived from available metadata. Use
% of metadata is recommended to retrieve site- & protocol-specific
% parameters and ensure appropriate data handling and processing.
% 
%
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium
%__________________________________________________________________________

global hmri_def

% =================== B1 mapping processing parameters ====================
% B1 map protocol parameters
% This is a template file that can be modified and saved with meaningful
% name in order to use customized B1 protocols 
% --------------------------

% 'i3D_AFI'
hmri_def.b1map.i3D_AFI.b1type = 'i3D_AFI'; 
hmri_def.b1map.i3D_AFI.b1avail = true; 
hmri_def.b1map.i3D_AFI.procreq = true; 
hmri_def.b1map.i3D_AFI.b1acq.TR2TR1ratio = 5;
hmri_def.b1map.i3D_AFI.b1acq.alphanom = 60;
hmri_def.b1map.i3D_AFI.deffnam  = 'hmri_b1_defaults.m';

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

% 'i3D_EPI'
hmri_def.b1map.i3D_EPI.b1type = 'i3D_EPI'; 
hmri_def.b1map.i3D_EPI.b1avail   = true; 
hmri_def.b1map.i3D_EPI.procreq = true; 
hmri_def.b1map.i3D_EPI.deffnam = 'hmri_b1_defaults.m';
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
hmri_def.b1map.i3D_EPI.b1acq.beta = 115:-5:65;
hmri_def.b1map.i3D_EPI.b1acq.TM = 31.2;
hmri_def.b1map.i3D_EPI.b1acq.tert = 540e-3*24; % EchoSpacing * numberPElines
hmri_def.b1map.i3D_EPI.b1acq.blipDIR = 1;
% b0-acquisition
hmri_def.b1map.i3D_EPI.b0acq.shortTE = 10; % ms
hmri_def.b1map.i3D_EPI.b0acq.longTE = 12.46; % ms
hmri_def.b1map.i3D_EPI.b0acq.iformat = 'PM'; % ms

% 'tfl_b1_map'
hmri_def.b1map.tfl_b1_map.b1type = 'tfl_b1_map'; 
hmri_def.b1map.tfl_b1_map.avail   = true; 
hmri_def.b1map.tfl_b1_map.procreq = true; 

% 'rf_map'
hmri_def.b1map.rf_map.b1type = 'rf_map'; 
hmri_def.b1map.rf_map.avail   = true; 
hmri_def.b1map.rf_map.procreq = true; 

end