% Sets the defaults for B1 bias correction, part of the hMRI toolbox.
% Consider this file as a template for local settings specifications. 
% Please read below for details.
%
% This local defaults file has settings which are useful for qMRI-BIDS
% compliant datasets.
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

% 'i3D_EPI'
% b1-validation
hmri_def.b1map.i3D_EPI.b1validation.checkTEs = true; % input validation using image TEs. Assumes SE has shorter TE than STE in metadata (qMRI-BIDS assumption). Disabled by default as this assumption is not valid for the metadata in the DICOMs from some sequences.
hmri_def.b1map.i3D_EPI.b1validation.useBidsFlipAngleField = true; % qMRI-BIDS support: read flip angles from BIDS metadata. Requires this metadata to be set appropriately!
