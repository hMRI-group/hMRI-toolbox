% Sets the defaults for B1 bias correction, part of the hMRI toolbox.
% Local settings to process the hMRI dataset available on
% https://github.molgen.mpg.de/hMRI-group/Toolbox.
%
% WARNING: Don't use these settings for other datasets without a careful
% check of the parameter values and their matching with your own data.
%
% FORMAT hmri_b1_local_defaults
%__________________________________________________________________________
%
% PURPOSE
% To set user-defined (site- or protocol-specific) defaults parameters for
% B1 mapping. Applies to 3D EPI, 3D AFI and UNICORT protocols only. 
% Customized processing parameters can be defined, overwriting defaults
% from hmri_b1_standard_defaults. Acquisition parameters can be specified
% here as a fallback solution when no metadata are available. Note that the
% use of metadata is strongly recommended.  
%
% RECOMMENDATIONS
% Parameters defined in this file are identical, initially, to the ones
% defined in hMRI-Toolbox\config\hmri_b1_standard_defaults. It is
% recommended, when modifying this file, to remove all unchanged entries
% and SAVE THE MODIFIED FILE WITH A MEANINGFUL NAME. This will help you
% identifying the appropriate defaults to be used during map creation for
% B1 map calculation, and will improve the readability of the file by
% pointing to the modified parameters only. 
%
% WARNING
% Modification of the defaults parameters may impair the integrity of the
% toolbox, leading to unexpected behaviour. Only recommended for expert
% users. 
%
% HOW DOES IT WORK?
% The modified defaults file can be selected when specifying the B1 type in
% the "Create maps" branch of the hMRI-Toolbox.
%__________________________________________________________________________
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium
%__________________________________________________________________________

% Global hmri_def variable used across the whole toolbox
global hmri_def

%--------------------------------------------------------------------------
% B1 mapping processing parameters 
%--------------------------------------------------------------------------

% 'i3D_EPI'
hmri_def.b1map.i3D_EPI.b1type = 'i3D_EPI'; 
hmri_def.b1map.i3D_EPI.b1avail   = true; 
hmri_def.b1map.i3D_EPI.procreq = true; 
% b0&b1-processing
hmri_def.b1map.i3D_EPI.b1proc.b0maskbrain = 1;
% b1-acquisition
hmri_def.b1map.i3D_EPI.b1acq.beta = 115:-5:65;
hmri_def.b1map.i3D_EPI.b1acq.TM = 33.8;
hmri_def.b1map.i3D_EPI.b1acq.tert = 540e-3*24; % EchoSpacing * numberPElines
hmri_def.b1map.i3D_EPI.b1acq.blipDIR = 1;
% b0-acquisition
hmri_def.b1map.i3D_EPI.b0acq.shortTE = 10; % ms
hmri_def.b1map.i3D_EPI.b0acq.longTE = 12.46; % ms
hmri_def.b1map.i3D_EPI.b0acq.iformat = 'PM'; % ms
