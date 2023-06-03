function hmri_b1_local_defaults_test_afi_bids
% Sets the defaults for B1 bias correction, part of the hMRI toolbox.
% Consider this file as a template for local settings specifications. 
% Please read below for details.
%
% This local defaults file disables smoothing for AFI in order to make results
% on synthetic data comparable with ground truth for testing.
%
% FORMAT hmri_b1_local_defaults

% Global hmri_def variable used across the whole toolbox
global hmri_def

% 'i3D_AFI'
hmri_def.b1map.i3D_AFI.b1proc.B1FWHM = 0; % For smoothing of B1 map. FWHM in mm; set to 0 to disable smoothing

end
