% These local settings have been defined for the hMRI toolbox for tutorial
% purpose only. Please read carefully the comments in the code below before
% reusing these parameters for your own processing.  

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

% The following examples can be un/commented for testing purpose.
% Parameters in examples 1-3 can be used and/or modified simultaneously or
% separately. 

%% EXAMPLE 1 
% To specify the name of the research centre where the data are
% processed. Although not mandatory, this parameter will be saved with all
% the defaults used for a given data processing. Therefore comes handy to
% trace the origin of a given set of generated maps.  
hmri_def.centre = 'my_research_centre' ; 

%% EXAMPLE 2
% To disable the cleanup of the intermediate files and directories. This
% comes handy to debug any processing bug leading to erroneous maps. Since
% all intermediate images are kept, the problem can more easily be tracked
% down to its origin. The cleanup option should be enabled otherwise, to
% limit the volume of data and number of images generated, for the sake of
% disk space and clarity.  
hmri_def.cleanup = false;

%% EXAMPLE 3
% To disable generation of QA results. When testing various options to
% process data, disabling the generation of QA results might be desirable
% to speed up the processing time. This option should not be disabled
% otherwise.
hmri_def.qMRI_maps.QA = 0; 

%% EXAMPLE 4
% Use (uncomment selectively) one of the following examples (a-c) to
% compare results obtained with different strategies for R2* bias
% correction (cf. Balteau et al., ISMRM 2018, p.2694).

% % EXAMPLE 4a - generate maps using OLS fit at TE=0 for each contrast
% % (default). 
% hmri_def.fullOLS = true;
% hmri_def.PDproc.T2scorr = 0;  

% % EXAMPLE 4b - generate maps using averaged echo of each contrast without
% % R2* bias correction (not recommended!).
% hmri_def.fullOLS = false;
% hmri_def.PDproc.T2scorr = 0;  

% % EXAMPLE 4c - generate maps using averaged echo of each contrast with R2*
% % bias correction.
% hmri_def.fullOLS = false;
% hmri_def.PDproc.T2scorr = 1;  
