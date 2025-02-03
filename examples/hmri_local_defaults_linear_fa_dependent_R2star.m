% PURPOSE
% To set user-defined (site- or protocol-specific) defaults parameters
% which are used by the hMRI toolbox. Customized processing parameters can
% be defined, overwriting defaults from hmri_defaults. Acquisition
% protocols can be specified here as a fallback solution when no metadata
% are available. Note that the use of metadata is strongly recommended. 
%
% RECOMMENDATIONS
% Parameters defined in this file are identical, initially, to the ones
% defined in hmri_defaults.m. It is recommended, when modifying this file,
% to remove all unchanged entries and save the file with a meaningful name.
% This will help you identifying the appropriate defaults to be used for
% each protocol, and will improve the readability of the file by pointing
% to the modified parameters only.
%
% WARNING
% Modification of the defaults parameters may impair the integrity of the
% toolbox, leading to unexpected behaviour. ONLY RECOMMENDED FOR ADVANCED
% USERS - i.e. who have a good knowledge of the underlying algorithms and
% implementation. The SAME SET OF DEFAULT PARAMETERS must be used to
% process uniformly all the data from a given study. 
%
% HOW DOES IT WORK?
% The modified defaults file can be selected using the "Configure toolbox"
% branch of the hMRI-Toolbox. For customization of B1 processing
% parameters, type "help hmri_b1_standard_defaults.m". 
%
% DOCUMENTATION
% A brief description of each parameter is provided together with
% guidelines and recommendations to modify these parameters. With few
% exceptions, parameters should ONLY be MODIFIED and customized BY ADVANCED
% USERS, having a good knowledge of the underlying algorithms and
% implementation. 
% Please refer to the documentation in the github WIKI for more details. 
%__________________________________________________________________________
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium

% Global hmri_def variable used across the whole toolbox
global hmri_def

%==========================================================================
% R1/PD/R2s/MT map creation parameters
%==========================================================================

% Choose model of flip-angle dependence of R2*
% - 'none': classic ESTATICS model assuming constant R2* between contrasts
% - 'linear': method from Melotti, et al. (MRM 2023). Note that this model
%             is not compatible with MT-weighted data.
hmri_def.R2s_flip_angle_dependence = 'linear';
